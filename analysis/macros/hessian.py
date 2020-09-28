#!/usr/bin/env python
from __future__ import print_function
import argparse
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

def _RooAbsCollection__iter__(self):
    it = self.iterator()
    obj = it.Next()
    while obj != None:
        yield obj
        obj = it.Next()

ROOT.RooAbsCollection.__iter__ = _RooAbsCollection__iter__

# because RooAbsCollection uses operator= for assignment :(
def rooAssign(target, other):
    if target == other:
        return
    for el in target:
        theirs = other.find(el)
        if not theirs:
            continue
        el.setVal(theirs.getVal())
        el.setError(theirs.getError())
        el.setAsymError(theirs.getErrorLo(), theirs.getErrorHi())
        el.setAttribute("Constant", theirs.isConstant())


def finitediff(f, x, y, dx, dy):
    x0 = x.getVal()
    y0 = y.getVal()

    # stencil
    def s(i, j):
        x.setVal(x0 + i*dx)
        y.setVal(y0 + j*dy)
        return f.getVal()

    out = 0.
    if x is y:
        # https://en.wikipedia.org/wiki/Five-point_stencil
        out += -s(2,2) + 16*s(1,1) - 30*s(0,0) + 16*s(-1,-1) - s(-2,-2)
        out /= 12*dx*dy
    else:
        # Following http://www.holoborodko.com/pavel/2014/11/04/computing-mixed-derivatives-by-finite-differences/
        out += 8*( s(1,-2) + s(2,-1) + s(-2,1) + s(-1,2) )
        out -= 8*( s(-1,-2) + s(-2,-1) + s(1,2) + s(2,1) )
        out -= ( s(2,-2) + s(-2,2) - s(-2,-2) - s(2,2) )
        out += 64*( s(-1,-1) + s(1,1) - s(1,-1) - s(-1,1) )
        out /= 144*dx*dy

    x.setVal(x0)
    y.setVal(y0)
    return out


def compute_hessian(args):
    fws, wsname = args.workspace.split(':')
    fin = ROOT.TFile.Open(fws)
    w = fin.Get(wsname)

    ffit, fitname = args.fit.split(':')
    fin2 = ROOT.TFile.Open(ffit)
    fit = fin2.Get(fitname)

    model = w.obj(args.model)
    pdf = model.GetPdf()
    nll = pdf.createNLL(w.data("data_obs"), ROOT.RooLinkedList())
    params = pdf.getParameters(model.GetObservables())
    rooAssign(params, fit.constPars())
    rooAssign(params, fit.floatParsFinal())

    floatparams = [p for p in params if fit.floatParsFinal().find(p)]
    npar = len(floatparams)
    hess = np.zeros(shape=(npar, npar))
    for ix in range(npar):
        x = floatparams[ix]
        dx = args.scale*x.getError()
        for iy in range(ix, npar):
            y = floatparams[iy]
            dy = args.scale*y.getError()
            hess[ix, iy] = finitediff(nll, x, y, dx, dy)

    ilow  = np.tril_indices(npar, -1)
    hess[ilow] = hess.T[ilow]
    return hess, floatparams, nll


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute hessian by hand.")
    parser.add_argument("-w", "--workspace", metavar="ROOTFILE:WORKSPACE", help="Workspace to load", required=True)
    parser.add_argument("-f", "--fit", metavar="ROOTFILE:FIT_NAME", help="Fit result to load", required=True)
    parser.add_argument("-m", "--model", help="Model to load", default="ModelConfig")
    parser.add_argument("-s", "--scale", help="Scale (multiplier of parameter error) to use for finite difference", default=0.5, type=float)
    parser.add_argument("--cond", help="Regularize the inversion by adding identity times this factor to the Hessian", default=0., type=float)

    args = parser.parse_args()
    hess, param, nll = compute_hessian(args)
    param = np.array(param)
    param_fit_err = np.fromiter((p.getError() for p in param), dtype='d')

    hess += np.eye(param.size) * args.cond
    cond = np.linalg.cond(hess)
    print("Condition number of hessian:", cond)
    hess_eig, hess_eigv = np.linalg.eigh(hess)

    normed_eigv = hess_eigv[0]
    order = np.abs(normed_eigv).argsort()[:-10:-1]
    print("Shallowest direction (eigenvalue %r, top 10):" % hess_eig[0])
    for mag, par in zip(normed_eigv[order], param[order]):
        print("  %6.3f %r" % (mag, par.GetName()))

    normed_eigv = hess_eigv[-1]
    order = np.abs(normed_eigv).argsort()[:-10:-1]
    print("Steepest direction (eigenvalue %r, top 10):" % hess_eig[-1])
    for mag, par in zip(normed_eigv[order], param[order]):
        print("  %6.3f %r" % (mag, par.GetName()))

    covar = np.linalg.inv(hess)
    cov_eig, cov_eigv = np.linalg.eigh(covar)
    print("Largest eigenvalue of covariance matrix:", cov_eig[-1])
    print("Smallest eigenvalue of covariance matrix:", cov_eig[0])

    std = np.sqrt(np.diag(covar))
    print("Largest standard deviation:", std.max(), "at", param[std.argmax()])
    print("Smallest standard deviation:", std.min(), "at", param[std.argmin()])

    corr = covar / std[:,None] / std[None,:]
    corr_eig, corr_eigv = np.linalg.eigh(corr)

    normed_eigv = corr_eigv[-1]
    order = np.abs(normed_eigv).argsort()[:-10:-1]
    print("Largest correlation vector (eigenvalue %r, top 10):" % corr_eig[-1])
    for mag, par in zip(normed_eigv[order], param[order]):
        print("  %6.3f %r" % (mag, par.GetName()))

