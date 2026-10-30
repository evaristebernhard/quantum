"""Quick invsqrt verification: full domain at 1 shot/input on 8 cores, plus top-input boundary check."""
from pathlib import Path
import os, sys, json, subprocess, time, math
RUN=Path(__file__).resolve().parent
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:8]))
ref=json.loads((RUN/'QCC_RELEASE_MANIFEST.json').read_text())['reference_costs']
env=os.environ.copy();env['PYTHONPATH']=str(RUN/'official_snapshot')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
env['PYTHONUNBUFFERED']='1'
# 3133440 is a multiple of batch 64; the remaining 12 highest inputs are checked separately below.
cmd=[sys.executable,'-m','invsqrt.verify_submission_invsqrt',str(RUN/'candidate/invsqrt.kmb'),
     '--shots-per-x','1','--max-x','3133440','--batch-size','64','--max-workers','8',
     '--salt','chayu-quick-20261030-a']
print('START_INV',json.dumps(cmd),flush=True)
t=time.monotonic()
try:
 p=subprocess.run(cmd,cwd=RUN,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=420)
 (RUN/'invsqrt-quick-aligned.stdout.txt').write_text(p.stdout)
 (RUN/'invsqrt-quick-aligned.stderr.txt').write_text(p.stderr)
 try:r=json.loads(p.stdout)
 except json.JSONDecodeError:r={'status':'PROCESS_ERROR','stdout_tail':p.stdout[-1500:],'stderr_tail':p.stderr[-1500:]}
 r['process_exit_code']=p.returncode
except subprocess.TimeoutExpired as ex:
 r={'status':'QUICK_TIMEOUT','limit_seconds':420}
 for label,data in [('stdout',ex.stdout),('stderr',ex.stderr)]:
  (RUN/f'invsqrt-quick-aligned.{label}.txt').write_bytes(data if isinstance(data,bytes) else (data or b''))
r['wall_seconds']=time.monotonic()-t
r['verification_mode']='1 shot per input x=1..3133440 plus 256-shot check of the 12 largest inputs; measurement outcomes sampled once each, not the official 64-shot full verification'
print('INV_DONE',json.dumps(r),flush=True)
if r.get('status')=='PASS':
 sys.path.insert(0,str(RUN/'official_snapshot'))
 import kickmix as km,numpy as np
 from invsqrt import verify_submission_invsqrt as iv
 c=km.Circuit.from_file(RUN/'candidate/invsqrt.kmb',format='kmb')
 spec=iv.ProblemSpecInvSqrt.for_n_p(10)
 inputs,scratch=iv.validate_registers(c,spec)
 xs=np.resize(np.arange(3133441,3133453,dtype=np.uint64),256)
 sim=km.Simulator(batch_size=256,seed=20261030,debug_print=False)
 sim.clear_for_shot();sim.write_across_shots(inputs,iv.to_bits(xs,spec.width));sim.do(c)
 assert np.array_equal(iv.from_bits(sim.read_across_shots(inputs)),xs),'x not restored on largest inputs'
 assert not sim.read_across_shots(scratch).any(),'scratch dirty on largest inputs'
 actual=np.array([iv.half_turns_to_units(sim.read_shot_phase(i)) for i in range(256)],dtype=np.uint64)
 roots=[math.isqrt((1<<130)//int(x)) for x in xs]
 target=np.array([((t+1)>>1)&((1<<64)-1) for t in roots],dtype=np.uint64)
 tol=np.array([min(t>>25,1<<63) for t in roots],dtype=np.uint64)
 err=np.abs((actual-target).view(np.int64)).view(np.uint64)
 assert (err<=tol).all(),'phase out of tolerance on largest inputs'
 r['highest_12_input_boundary_check']={'status':'PASS','shots':256,'worst_error_over_tolerance':float((err/tol).max())}
 r['quick_cost_ratio']=ref['invsqrt']/r['score']
(RUN/'invsqrt-quick-aligned.json').write_text(json.dumps(r,indent=2))
# Consolidated quick summary with an explicit, clearly-labeled projection for DQI.
summary=json.loads((RUN/'quick-summary.json').read_text())
if r.get('status')=='PASS':summary['reports']['invsqrt']=r
summary['quick_gm_3_tasks']=None
reports=summary['reports']
if all(reports.get(k,{}).get('status')=='PASS' for k in ['invsqrt','dqi','dfthcsos']):
 summary['quick_gm_3_tasks']=math.prod(reports[k]['quick_cost_ratio'] for k in reports)**(1/3)
 from dqi import verify_submission_dqi as dv
 d=reports['dqi'];eps=dv.clopper_pearson_upper_bound(0,102400,log2_delta=-128)
 projected_cost=dv.repetition_overhead(eps)*d['num_qubits']**2*d['toffoli_count']
 summary['dqi_full_shot_projection']={'assumption':'ZERO errors in a hypothetical 102400-shot run at the same mean Toffoli count; NOT verified and not a leaderboard score.','dqi_cost':projected_cost,'dqi_ratio':ref['dqi']/projected_cost}
 if reports['invsqrt'].get('quick_cost_ratio'):
  summary['projected_gm_if_dqi_zero_error']=(reports['invsqrt']['quick_cost_ratio']*reports['dfthcsos']['quick_cost_ratio']*(ref['dqi']/projected_cost))**(1/3)
summary['full_verification']=False
(RUN/'quick-summary-final.json').write_text(json.dumps(summary,indent=2))
print('FINAL_SUMMARY',json.dumps(summary),flush=True)
