"""Run unmodified official verifier CLIs with explicitly reduced shots."""
from pathlib import Path
import os, sys, json, subprocess, time, math, hashlib
RUN=Path(__file__).resolve().parent
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
start=time.monotonic()
ref=json.loads((RUN/'QCC_RELEASE_MANIFEST.json').read_text())['reference_costs']
reports={}
plan=[
 ('invsqrt',['-m','invsqrt.verify_submission_invsqrt',str(RUN/'candidate/invsqrt.kmb'),'--shots-per-x','1','--max-workers','4','--salt','chayu-quick-20261030-a'],1,360),
 ('dfthcsos',['-m','dfthcsos.verify',str(RUN/'candidate/dfthcsos'),'--shots-per-address','1','--salt','chayu-quick-20261030-a'],4,360),
 ('dqi',['-m','dqi.verify_submission_dqi',str(RUN/'candidate/dqi.kmb'),'--num-shots','4096','--max-workers','4','--salt','chayu-quick-20261030-a'],1,360),
]
for name,args,threads,limit in plan:
 env=os.environ.copy();env['PYTHONPATH']=str(RUN/'official_snapshot')
 for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]=str(threads)
 env['PYTHONUNBUFFERED']='1'
 cmd=[sys.executable]+args
 print('START',name, json.dumps(cmd),flush=True)
 t=time.monotonic()
 try:
  p=subprocess.run(cmd,cwd=RUN,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=limit)
  (RUN/f'{name}-quick.stdout.txt').write_text(p.stdout)
  (RUN/f'{name}-quick.stderr.txt').write_text(p.stderr)
  try:report=json.loads(p.stdout)
  except json.JSONDecodeError:report={'status':'PROCESS_ERROR','stdout_tail':p.stdout[-2000:],'stderr_tail':p.stderr[-2000:]}
  report['process_exit_code']=p.returncode
 except subprocess.TimeoutExpired as ex:
  report={'status':'QUICK_TIMEOUT','limit_seconds':limit}
  for label,data in [('stdout',ex.stdout),('stderr',ex.stderr)]:
   (RUN/f'{name}-quick.{label}.txt').write_bytes(data if isinstance(data,bytes) else (data or '').encode())
 report['wall_seconds']=time.monotonic()-t
 report['verification_mode']='reduced-shot quick check, not full official verification'
 if report.get('status')=='PASS':report['quick_cost_ratio']=ref[name]/report['score']
 (RUN/f'{name}-quick.json').write_text(json.dumps(report,indent=2))
 reports[name]=report
 summary={'verification_mode':'QUICK_ONLY','cpu_affinity':sorted(os.sched_getaffinity(0)),'total_seconds':time.monotonic()-start,'reports':reports,'warning':'DQI confidence and repetition overhead depend on sample count; quick GM is NOT the leaderboard metric or a certified prediction.'}
 if len(reports)==3 and all(r.get('status')=='PASS' for r in reports.values()):
  summary['reduced_shot_gm'] = math.prod(r['quick_cost_ratio'] for r in reports.values())**(1/3)
 (RUN/'quick-summary.json').write_text(json.dumps(summary,indent=2))
 print('DONE',name,json.dumps(report),flush=True)
print('QUICK_SUMMARY',json.dumps(summary),flush=True)
