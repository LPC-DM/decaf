selections_Run2 = dict()

selections_Run2['2016'] = dict()
selections_Run2['2017'] = dict()
selections_Run2['2018'] = dict()

### Basic 2016 selections

selections_Run2['2016']['ele_loose_id'] = 'Electron_mvaSpring16GP_WP90'
selections_Run2['2016']['ele_tight_id'] = 'Electron_mvaSpring16GP_WP80'
selections_Run2['2016']['ele_loose_pt'] = 7.0
selections_Run2['2016']['ele_tight_pt'] = 30.0
selections_Run2['2016']['ele_loose_eta'] = 2.4
selections_Run2['2016']['ele_loose_dxy'] = 0.05
selections_Run2['2016']['ele_loose_dz'] = 0.2
selections_Run2['2016']['ele_loose_iso'] = 0.4
selections_Run2['2016']['ele_tight_iso'] = 0.06
selections_Run2['2016']['ele_isovar'] = 'Electron_pfRelIso03_all'

selections_Run2['2016']['muo_loose_pt'] = 5.0
#selections_Run2['2016']['muo_tight_pt'] = 30.0
selections_Run2['2016']['muo_loose_eta'] = 2.4
selections_Run2['2016']['muo_loose_dxy'] = 0.5
selections_Run2['2016']['muo_loose_dz'] = 1.0
selections_Run2['2016']['muo_loose_iso'] = 0.4

selections_Run2['2016']['tau_loose_pt'] = 18.0
#selections_Run2['2016']['tau_tight_pt'] = 30.0
selections_Run2['2016']['tau_loose_eta'] = 2.3
selections_Run2['2016']['tau_loose_id'] = 'Tau_idMVAnew'

selections_Run2['2016']['pho_loose_pt'] = 15.0
#selections_Run2['2016']['pho_tight_pt'] = 30.0
selections_Run2['2016']['pho_loose_eta'] = 2.5

selections_Run2['2016']['fjet_pt'] = 200.0
selections_Run2['2016']['fjet_eta'] = 2.4
selections_Run2['2016']['fjet_deltaR'] = 1.5

selections_Run2['2016']['jet_pt'] = 25.0
selections_Run2['2016']['jet_eta'] = 4.5
selections_Run2['2016']['jet_deltaR'] = 0.4

selections_Run2['2016']['recoil_pt'] = 200.0
selections_Run2['2016']['recoil_jetpt_skinny'] = 100.0
selections_Run2['2016']['recoil_jetpt_loose'] = 200.0
selections_Run2['2016']['recoil_minDPhi_skinny'] = 0.5
selections_Run2['2016']['recoil_minDPhi_loose'] = 0.8

### Copy 2016 to 2017 and 2018
selections_Run2['2017'] = selections_Run2['2016'].copy()
selections_Run2['2018'] = selections_Run2['2016'].copy()

### Make 2017 changes
selections_Run2['2017']['ele_loose_id'] = 'Electron_mvaFall17Iso_WP90'
selections_Run2['2017']['ele_tight_id'] = 'Electron_mvaFall17Iso_WP80'
selections_Run2['2017']['ele_loose_iso'] = float('Infinity') # Isolation is baked into the ID in 2017
selections_Run2['2017']['ele_tight_iso'] = float('Infinity')
#selections_Run2['2017']['tau_loose_id'] = 'Tau_idMVAoldDM2017v2' # For 2017 NanoAODv2
selections_Run2['2017']['tau_loose_id'] = 'Tau_idAntiEle' # FIXME: should have a complete ID
