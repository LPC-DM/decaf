grep -i -l "error" logs/Z_boson_pt_* | sed s#logs#../../useful_scripts/scripts#g | sed 's/sh.*/sh/'
grep -i -l "error" logs/Z_boson_pt_* | sed s#logs#../../useful_scripts/scripts#g | sed 's/sh.*/sh/' | tr '\n' ' '
grep -i -l "error" logs/Z_boson_pt_*.err | sed s#logs#../../useful_scripts/scripts#g | sed 's/_26.*err/.sh/g' | tr '\n' ' '
for i in {0..1637};
    if [ ! -f Z_boson_pt_$i.root ]; then
        echo $i
    fi
done;


grep -i -L "finished" logs/Zll_2018_submitScript.13567097_*.out | xargs grep -h "sh" | xargs python NAFSubmit.py -r 360 -M 1000 -u -v /nfs/dust/cms/user/mwassmer/proxy/x509up_u26621
for i in $(grep -i -L "finished" logs/Zll_2018_submitScript.13567097_*.out | xargs); do ls ${${i}%.out}.*;done
