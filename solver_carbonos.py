import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

N=168 #total variables for a 24 hour plan

CORE=['steam_demand_tph','waste_heat_available_tph','fixed_electricity_mw',
      'baseline_flexible_mw','solar_available_mw','grid_factor_tco2e_per_mwh','tariff_inr_per_mwh']

def read_inputs(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames):
            raise ValueError('Missing or duplicate CSV headers')
        missing=set(['hour']+CORE)-set(reader.fieldnames)
        if missing:raise ValueError('Missing CSV columns: '+', '.join(sorted(missing)))
        rows=list(reader)
    if len(rows)!=24:raise ValueError('CSV must contain exactly 24 rows')
    parsed={}
    for row in rows:
        h=float(row['hour'])
        if not np.isfinite(h) or not h.is_integer() or not 0<=h<24:raise ValueError('hour must be an integer from 0 to 23')
        h=int(h)
        if h in parsed:raise ValueError('Duplicate hour: '+str(h))
        r={key:float(row[key]) for key in CORE}
        defaults={'flexible_min_mw':.5 if 8<=h<18 else 0,
                  'flexible_max_mw':2.5 if 6<=h<23 else 0,
                  'biomass_max_tph':65,'grid_max_mw':22,'duration_h':1}
        for key,default in defaults.items():r[key]=float(row[key]) if key in row else default
        if not all(np.isfinite(v) and v>=0 for v in r.values()):raise ValueError('Nonfinite/negative input at hour '+str(h))
        if r['duration_h']!=1:raise ValueError('Only one-hour intervals are supported')
        if not 20<=r['biomass_max_tph']<=65:raise ValueError('Biomass max must be 20 to 65 for the continuously-on model')
        if not 0<=r['grid_max_mw']<=22:raise ValueError('Grid limit exceeds case bounds')
        lo,hi=r['flexible_min_mw'],r['flexible_max_mw']
        if not 0<=lo<=hi<=defaults['flexible_max_mw'] or lo<defaults['flexible_min_mw']:
            raise ValueError('Flexible bounds violate production rules at hour '+str(h))
        parsed[h]=r
    return {key:np.array([parsed[h][key] for h in range(24)]) for key in parsed[0]}

def activate(data):
    global D,W,F,BASEF,S,EF,T,FMIN,FMAX,BMAX,GRIDMAX,C,E
    D,W,F,BASEF,S,EF,T=[data[k].copy() for k in CORE]
    FMIN=data['flexible_min_mw'];FMAX=data['flexible_max_mw']
    BMAX=data['biomass_max_tph'];GRIDMAX=data['grid_max_mw']
    C=np.concatenate([np.full(24,1530),np.full(24,1485),np.full(24,2560),np.zeros(48),T,np.zeros(24)])
    E=np.concatenate([np.full(24,.423),np.full(24,.027),np.full(24,.16),np.zeros(48),EF,np.zeros(24)])

def solve(carbon_cap=None,cost_cap=None,objective='cost'): #with no distrubances
    solar=S.copy(); bm=BMAX.copy()
    bounds=[]
    for lo,hi in [(25,np.full(24,90)),(20,bm),(0,np.full(24,45)),(0,W),(0,solar),(0,GRIDMAX)]:
        bounds.extend([(lo,float(v)) for v in hi])
    bounds.extend([(float(FMIN[h]),float(FMAX[h])) for h in range(24)])
    A=[];b=[];Ae=[];be=[]
    def row(items):
        v=np.zeros(N)
        for i,a in items:v[i]=a
        return v
    for h in range(24):
        r=row([(24*j+h,1) for j in range(4)])
        A.extend([r,-r]);b.extend([1.03*D[h],-D[h]])
        Ae.append(row([(96+h,1),(120+h,1),(144+h,-1)]));be.append(F[h])
    Ae.append(row([(144+h,1) for h in range(24)]));be.append(18)
    for j,ramp,initial in [(0,10,60),(1,8,20),(2,45,0),(3,20,8)]:
        for h in range(24):
            r=row([(24*j+h,1)]+([(24*j+h-1,-1)] if h else []))
            A.extend([r,-r]);b.extend([ramp+(initial if h==0 else 0),ramp-(initial if h==0 else 0)])
    for j,coef,limit in [(0,.18,350),(1,.27,220),(2,80,20000)]:
        A.append(row([(24*j+h,coef) for h in range(24)]));b.append(limit)
    if carbon_cap is not None:A.append(E);b.append(carbon_cap)
    if cost_cap is not None:A.append(C);b.append(cost_cap)
    if objective=='smooth':
        A=[np.r_[v,np.zeros(48)] for v in A];Ae=[np.r_[v,np.zeros(48)] for v in Ae]
        for j,initial in [(0,60),(1,20)]:
            for h in range(24):
                v=np.zeros(N+48);v[j*24+h]=1
                if h:v[j*24+h-1]=-1
                for sign in [1,-1]:
                    r=sign*v.copy();r[N+j*24+h]=-1;A.append(r);b.append(sign*initial if h==0 else 0)
        out=linprog(np.r_[np.zeros(N),np.ones(48)],A_ub=A,b_ub=b,A_eq=Ae,b_eq=be,bounds=bounds+[(0,None)]*48,method='highs')
        if out.success:out.x=out.x[:N]
    else:
        out=linprog(C if objective=='cost' else E,A_ub=A,b_ub=b,A_eq=Ae,b_eq=be,bounds=bounds,method='highs')
    if not out.success:raise RuntimeError('No optimal solution: '+out.message)
    return out

def audit(x):
    q,b,g,w,s,i,f=x.reshape(7,24);sa=S.copy();bmax=BMAX.copy()
    checks={
      'steam_min':np.min(q+b+g+w-D),'steam_max':np.min(1.03*D-q-b-g-w),
      'electric_balance_error':-np.max(np.abs(i+s-F-f)),
      'coal_bounds':min(np.min(q-25),np.min(90-q)),
      'bio_bounds':min(np.min(b-20),np.min(bmax-b)),
      'gas_bounds':min(np.min(g),np.min(45-g)),
      'waste_bounds':min(np.min(w),np.min(W-w)),
      'solar_bounds':min(np.min(s),np.min(sa-s)),
      'grid_bounds':min(np.min(i),np.min(GRIDMAX-i)),
      'flex_bounds':min(np.min(f-FMIN),np.min(FMAX-f)),
      'flex_forbidden':-np.max(np.abs(f[[0,1,2,3,4,5,23]])),
      'flex_daily_error':-abs(sum(f)-18),'continuity':np.min(f[8:18]-.5),
      'coal_fuel_margin':350-.18*sum(q),'bio_fuel_margin':220-.27*sum(b),
      'gas_fuel_margin':20000-80*sum(g),
    }
    for name,v,initial,ramp in [('coal',q,60,10),('bio',b,20,8),('gas',g,0,45),('waste',w,8,20)]:
        checks[name+'_ramp_margin']=ramp-np.max(np.abs(np.diff(np.r_[initial,v])))
    if min(checks.values()) < -1e-6:
        raise ValueError('Independent feasibility audit failed: '+str(checks))
    return {k:float(v) for k,v in checks.items()}

def summary(x):
    q,b,g,w,s,i,f=x.reshape(7,24)
    return dict(cost=float(C@x),emissions=float(E@x),cost_saving_pct=float(100*(1-C@x/5526720)),carbon_saving_pct=float(100*(1-E@x/986.992)),coal_t=float(.18*sum(q)),biomass_t=float(.27*sum(b)),gas_sm3=float(80*sum(g)),steam_totals=[float(sum(v)) for v in [q,b,g,w]],solar_mwh=float(sum(s)),grid_mwh=float(sum(i)),flex_mwh=float(sum(f)))


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)

def check_baseline():
    wh=np.minimum(W,10);co=np.minimum(90,D-wh-20);ga=D-wh-co-20
    x=np.concatenate([co,np.full(24,20),ga,wh,.8*S,F+BASEF-.8*S,BASEF])
    names=['steam_t','coal_steam_t','biomass_steam_t','gas_steam_t','waste_steam_t',
           'grid_mwh','solar_mwh','steam_emissions','grid_emissions','total_emissions',
           'steam_cost','grid_cost','total_cost']
    actual=[sum(D),sum(co),480,sum(ga),sum(wh),sum(x[120:144]),sum(.8*S),
            E[:96]@x[:96],EF@x[120:144],E@x,C[:96]@x[:96],T@x[120:144],C@x]
    official=[2797,1886,480,205,226,216.5,59.2,843.538,143.454,986.992,4123180,1403540,5526720]
    rows=[dict(metric=k,calculated=float(a),official=o,difference=float(a-o),
               pass_within_1pct=bool(abs(a-o)<=.01*o)) for k,a,o in zip(names,actual,official)]
    if min(x)<-1e-6 or not all(r['pass_within_1pct'] for r in rows):
        raise ValueError('Baseline differs from Table 6 by more than 1%; check normal CSV data')
    return rows

def main():
    p=argparse.ArgumentParser(description='Optimise one input CSV without changing its data.')
    p.add_argument('csv_file',type=Path)
    p.add_argument('--carbon-limit',type=float,help='Optional daily tCO2e limit; normal qualification uses 868.55296')
    p.add_argument('--check-baseline',action='store_true',help='Verify Table 6 from undisturbed inputs')
    p.add_argument('--out',type=Path,help='Output folder; default is <CSV filename>_results')
    args=p.parse_args()
    if args.carbon_limit is not None and (not np.isfinite(args.carbon_limit) or args.carbon_limit<0):
        raise ValueError('Carbon limit must be finite and nonnegative')
    activate(read_inputs(args.csv_file))
    baseline_rows=check_baseline() if args.check_baseline else None
    first=solve(carbon_cap=args.carbon_limit)
    second=solve(carbon_cap=args.carbon_limit,cost_cap=float(C@first.x)+1e-7,objective='carbon')
    limit=float(E@second.x)+1e-7
    if args.carbon_limit is not None:limit=min(limit,args.carbon_limit)
    third=solve(carbon_cap=limit,cost_cap=float(C@first.x)+1e-7,objective='smooth')
    x=third.x.copy();x[np.abs(x)<1e-6]=0
    checks=audit(x);result=summary(x)
    if abs(result['cost']-float(C@first.x))>1e-3 or result['emissions']>limit+1e-6:
        raise ValueError('Final schedule failed objective/limit verification')
    result.update(steam_emissions=float(E[:96]@x[:96]),grid_emissions=float(EF@x[120:144]),
                  steam_cost=float(C[:96]@x[:96]),grid_cost=float(T@x[120:144]),
                  carbon_gate_pass=bool(E@x<=868.55296+1e-6),cost_gate_pass=bool(C@x<=5250384+1e-6))
    rows=[];v=x.reshape(7,24)
    for h in range(24):
        q,b,g,w,s,i,f=map(float,v[:,h])
        margins={'steam_min':q+b+g+w-D[h],'steam_max':1.03*D[h]-q-b-g-w,
                 'coal_min':q-25,'coal_max':90-q,'bio_min':b-20,'bio_max':BMAX[h]-b,
                 'gas_min':g,'gas_max':45-g,'waste_min':w,'waste_available':W[h]-w,
                 'solar_min':s,'solar_available':S[h]-s,'grid_min':i,'grid_max':GRIDMAX[h]-i,
                 'flex_min':f-FMIN[h],'flex_max':FMAX[h]-f}
        for j,name,initial,ramp in [(0,'coal',60,10),(1,'bio',20,8),(2,'gas',0,45),(3,'waste',8,20)]:
            margins[name+'_ramp']=float(ramp-abs(v[j,h]-(v[j,h-1] if h else initial)))
        row=dict(hour=h,coal_steam_tph=q,biomass_steam_tph=b,gas_steam_tph=g,waste_heat_steam_tph=w,
                 solar_used_mw=s,grid_import_mw=i,flexible_load_mw=f,
                 emissions_tco2e=float(E.reshape(7,24)[:,h]@v[:,h]),cost_inr=float(C.reshape(7,24)[:,h]@v[:,h]),
                 electricity_balance_error_mw=float(i+s-F[h]-f),constraint_status='PASS',
                 binding_constraints='; '.join(k for k,val in margins.items() if abs(val)<1e-5))
        row.update({k:float(val) for k,val in margins.items()});rows.append(row)
    report={'input_csv':str(args.csv_file),'carbon_limit':args.carbon_limit,'summary':result,
            'assumptions':'Coal and biomass continuously ON; minimum-on requirements satisfied; starts=0; 1h intervals.',
            'note':'Gate flags compare to the official normal thresholds; they are not required for the disturbance case.',
            'solver_stages':[r.message for r in [first,second,third]],'audit_tolerance':1e-6,
            'daily_constraint_checks':checks,'baseline_verification':baseline_rows}
    folder=args.out or Path(args.csv_file.stem+'_results');folder.mkdir(parents=True,exist_ok=True)
    write_csv(folder/'hourly_schedule.csv',rows)
    if baseline_rows:write_csv(folder/'baseline_verification.csv',baseline_rows)
    (folder/'results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    (folder/'report.md').write_text('# CarbonOS results\n\n'+report['assumptions']+'\n\n'+report['note']+
        '\n\n```json\n'+json.dumps(result,indent=2)+'\n```\n\n## Feasibility checks\n\n```json\n'+
        json.dumps(checks,indent=2)+'\n```\n\nFull hourly set-points, costs, emissions and binding constraints: hourly_schedule.csv.\n',encoding='utf-8')
    print(f"Cost: INR {result['cost']:,.2f} | Saving: {result['cost_saving_pct']:.4f}%")
    print(f"Emissions: {result['emissions']:.6f} tCO2e | Reduction: {result['carbon_saving_pct']:.4f}%")
    print('Feasibility: PASS. Reports:',folder.resolve())

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,OSError,KeyError,TypeError) as error:
        print('Error:',error,file=sys.stderr);sys.exit(1)