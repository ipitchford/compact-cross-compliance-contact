"""Reject plausible damaged algebra in isolated modules, also under Python -O."""
import pathlib
import subprocess
import sys
import tempfile

root=pathlib.Path(__file__).resolve().parent
cases=[('verify.py','age_peak_lower = F(3, 2) + F(4, 3)','age_peak_lower = F(1, 2)','exact_sign_certificates'),
       ('revision_checks.py','cm = sp.Rational(1,2)/(1+2)','cm = sp.Rational(1,2)/(1+1)','exact_elastic_checks'),
       ('revision_checks.py','delta == -sp.Rational(7,40)','delta == sp.Rational(7,40)','exact_elastic_checks')]
for filename,old,new,call in cases:
    original=(root/filename).read_text()
    if original.count(old)!=1:
        raise RuntimeError('Mutation anchor not unique: '+old)
    with tempfile.TemporaryDirectory(prefix='friction-negative-') as tmp:
        p=pathlib.Path(tmp)/'mutant.py';p.write_text(original.replace(old,new))
        for flags in ([],['-O']):
            code='import runpy; m=runpy.run_path('+repr(str(p))+'); m['+repr(call)+']()'
            done=subprocess.run([sys.executable,*flags,'-c',code],capture_output=True,text=True)
            if done.returncode==0 or 'AssertionError' not in done.stderr:
                raise RuntimeError('Corruption was not rejected: '+old)
            print('REJECT',filename,call,'optimized='+str(bool(flags)))
print('PASS: six semantic mutation invocations rejected')
