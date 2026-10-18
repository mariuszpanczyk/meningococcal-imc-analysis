"""IMD/IMC code-only companion, v0.2.0. Public mode does not load study data.
See the notebook and README for scientific provenance and limitations.
"""

INSTALL_REFERENCE_ENVIRONMENT = False
if INSTALL_REFERENCE_ENVIRONMENT:
    import subprocess, sys
    subprocess.check_call([sys.executable,'-m','pip','install','--quiet',
        'numpy==2.1.3','pandas==2.2.3','scipy==1.14.1',
        'statsmodels==0.14.4','matplotlib==3.9.2','nbformat==5.10.4'])

import io, json, hashlib, warnings, platform, re, unicodedata
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np
import pandas as pd
from scipy import stats, optimize, special
from scipy.stats.contingency import odds_ratio
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RUN_PRIVATE_ANALYSIS = False
DATA_PATH = None
ALLOW_COLAB_UPLOAD = False
CONFIRM_AUTHORIZED_ENVIRONMENT = False
CONFIRM_ANALYSIS_COHORT_ELIGIBILITY = False
FULL_COHORT_PATH = None  # Optional separate file; never infer missing participants.
EXPORT_PRIVATE_AGGREGATES = False
PRIVATE_OUTPUT_DIR = Path('private_results')
SEED = 20261018  # New reconstruction seed, not a recovered historical seed.
N_PERM = 19999
N_BOOT = 4999
RUN_FIRTH = True
RUN_SUPPLEMENTARY_WY = True

Y = 'powiklanie_immuno'
CLINICAL_NUM = ['wiek_lata','objawy_do_hosp_dni','przyjecie_do_atb_dni']
CLINICAL_CAT = ['plec','zomr','serotyp']
LAB = ['crp_max','pct_max','wbc_max','wbc_min','limf_min','plt_min',
       'mleczany_max','hco3_min','inr_max','fibrynogen_min','d_dimer_max']
CATEGORICAL = {'plec':['K','M'],'zomr':['nie','nie badano','tak'],
               'serotyp':['B','C','W','Y'],'oitd':['nie','tak']}
FAMILIES = {'clinical':CLINICAL_NUM+CLINICAL_CAT,'laboratory_severity':LAB+['oitd']}
LOG2 = {'crp_max','pct_max','mleczany_max','d_dimer_max'}
DATE_COLUMNS = ['data_urodzenia','data_objawow','data_przyjecia','data_atb','data_wypisu']
OPTIONAL_NUM = ['czas_hosp_dni','dlugosc_atb_dni']
UNITS = {'wiek_lata':'years','objawy_do_hosp_dni':'days','przyjecie_do_atb_dni':'days',
 'crp_max':'mg/dL','pct_max':'ng/mL','wbc_max':'10^3/uL','wbc_min':'10^3/uL',
 'limf_min':'10^3/uL','plt_min':'10^3/uL','mleczany_max':'mmol/L',
 'hco3_min':'mmol/L','inr_max':'unitless','fibrynogen_min':'g/L','d_dimer_max':'ng/mL'}

def rng_for(label):
    k=int.from_bytes(hashlib.sha256(label.encode()).digest()[:8],'little')
    return np.random.default_rng(np.random.SeedSequence([SEED,k]))

def bh(p):
    a=np.asarray(p,float); out=np.full(a.shape,np.nan); ok=np.isfinite(a)
    if ok.any(): out[ok]=multipletests(a[ok],method='fdr_bh')[1]
    return out

def clean_text(x):
    return unicodedata.normalize('NFKC',str(x)).strip()

def number(s):
    t=s.astype('string').str.strip().str.replace(',','.',regex=False)
    t=t.mask(t.str.lower().isin(['','na','nan','n/a','none','null']))
    try: z=pd.to_numeric(t,errors='raise').astype(float)
    except Exception: raise ValueError('Invalid numeric input; inspect the private source.') from None
    if np.isinf(z.to_numpy()).any(): raise ValueError('Infinite numeric input.')
    return z

def read_private(path=None,allow_upload=False):
    if path is None:
        if not allow_upload: raise ValueError('Set DATA_PATH or explicitly authorize a Colab upload.')
        from google.colab import files
        uploaded=files.upload()
        if len(uploaded)!=1: raise ValueError('Upload exactly one CSV file.')
        payload=next(iter(uploaded.values()))
    else: payload=Path(path).read_bytes()
    raw=pd.read_csv(io.BytesIO(payload),sep=';',encoding='utf-8-sig',dtype=str,skip_blank_lines=True)
    raw.columns=[clean_text(c) for c in raw.columns]
    if raw.columns.duplicated().any(): raise ValueError('Duplicate column names.')
    raw=raw.replace(r'^\s*$',np.nan,regex=True).dropna(how='all').reset_index(drop=True)
    return raw,hashlib.sha256(payload).hexdigest()

def prepare(raw,full_cohort=False):
    required=['plec','wiek_lata','objawy_do_hosp_dni','zomr','oitd','serotyp',Y]+LAB
    absent=[c for c in required if c not in raw]
    if absent: raise ValueError('Missing required columns: '+', '.join(absent))
    audit=[]
    if 'id' in raw:
        if raw['id'].isna().any() or raw['id'].duplicated().any():
            raise ValueError('Missing or duplicate private participant identifiers.')
        audit.append({'check':'duplicate_ID','violations':0})
    else: audit.append({'check':'duplicate_ID_not_assessed','violations':np.nan})
    d=pd.DataFrame(index=raw.index)
    for c in ['wiek_lata','objawy_do_hosp_dni']+LAB+[c for c in OPTIONAL_NUM if c in raw]:
        d[c]=number(raw[c])
        if d[c].dropna().lt(0).any(): raise ValueError('Negative numeric measurement: '+c)
    for c,levels in CATEGORICAL.items():
        s=raw[c].map(lambda v: clean_text(v) if pd.notna(v) else np.nan)
        if c in ['zomr','oitd']: s=s.str.lower()
        else: s=s.str.upper()
        if not s.dropna().isin(levels).all(): raise ValueError('Unknown category in '+c)
        d[c]=s
    outcome=raw[Y].str.strip().str.lower()
    if not outcome.dropna().isin(['tak','nie']).all(): raise ValueError('Unknown outcome code.')
    d[Y]=outcome.map({'nie':0.,'tak':1.})
    if not full_cohort and d[Y].isna().any(): raise ValueError('Analysis cohort requires observed IMC status.')
    if not full_cohort and set(d[Y].unique())!={0.,1.}: raise ValueError('Both outcome groups required.')
    dates={}
    for c in DATE_COLUMNS:
        if c in raw:
            dates[c]=pd.to_datetime(raw[c],format='%d.%m.%Y',errors='coerce')
            if (raw[c].notna() & dates[c].isna()).any(): raise ValueError('Invalid date format in '+c)
    if 'data_atb' in dates and 'data_przyjecia' in dates:
        delay=(dates['data_atb']-dates['data_przyjecia']).dt.total_seconds()/86400
        if 'przyjecie_do_atb_dni' in raw:
            provided=number(raw['przyjecie_do_atb_dni'])
            both=provided.notna() & delay.notna()
            if not np.allclose(provided[both],delay[both]): raise ValueError('Conflicting antibiotic delay.')
        d['przyjecie_do_atb_dni']=delay
    elif 'przyjecie_do_atb_dni' in raw: d['przyjecie_do_atb_dni']=number(raw['przyjecie_do_atb_dni'])
    else: raise ValueError('Provide admission/antibiotic dates or explicit admission-to-antibiotic delay.')
    if d['przyjecie_do_atb_dni'].dropna().lt(0).any(): raise ValueError('Negative admission-to-antibiotic delay.')
    for a,b in [('data_urodzenia','data_objawow'),('data_objawow','data_przyjecia'),
                ('data_przyjecia','data_atb'),('data_atb','data_wypisu')]:
        if a in dates and b in dates:
            count=int((dates[a]>dates[b]).sum())
            audit.append({'check':a+'_before_'+b,'violations':count})
            if count: raise ValueError('Chronology violation; inspect source privately.')
    if 'data_objawow' in dates and 'data_przyjecia' in dates:
        calc=(dates['data_przyjecia']-dates['data_objawow']).dt.days
        both=calc.notna() & d['objawy_do_hosp_dni'].notna()
        if not np.allclose(calc[both],d.loc[both,'objawy_do_hosp_dni']):
            raise ValueError('Conflicting onset-to-admission delay.')
    if 'oczekiwanie_atb_dni' in raw and 'data_objawow' in dates and 'data_atb' in dates:
        calc=(dates['data_atb']-dates['data_objawow']).dt.days
        supplied=number(raw['oczekiwanie_atb_dni']); both=calc.notna() & supplied.notna()
        audit.append({'check':'source_waiting_time_vs_onset_to_antibiotic',
                      'violations':int((~np.isclose(calc[both],supplied[both])).sum())})
    for c in LAB:
        x=d[c].dropna(); q1,q3=x.quantile([.25,.75]); iqr=q3-q1
        audit.append({'check':'IQR_flags_'+c,'violations':int(((x<q1-1.5*iqr)|(x>q3+1.5*iqr)).sum())})
    if full_cohort:
        if 'death_within_24h' not in raw or 'data_objawow' not in dates:
            raise ValueError('Full-cohort annual plot requires death_within_24h and onset dates.')
        death=raw['death_within_24h'].str.strip().str.lower().map({'tak':1.,'nie':0.})
        if death.isna().any(): raise ValueError('Explicit yes/no early-death status required.')
        if (d[Y].isna() & death.ne(1)).any(): raise ValueError('Unassessable IMC not explained by early-death status.')
        if (death.eq(1) & d[Y].notna()).any(): raise ValueError('Early-death/IMC eligibility conflict.')
        d['death_within_24h']=death
        d['onset_year']=dates['data_objawow'].dt.year
        if d['onset_year'].isna().any(): raise ValueError('Missing onset year.')
    return d,pd.DataFrame(audit)

def descriptive(d):
    rows=[]
    for c in CLINICAL_NUM+LAB+[v for v in OPTIONAL_NUM if v in d]:
        for g in [0,1]:
            x=d.loc[d[Y].eq(g),c].dropna()
            rows.append({'variable':c,'IMC':g,'n':len(x),'median':x.median(),
                         'Q1':x.quantile(.25),'Q3':x.quantile(.75),'unit':UNITS.get(c,'days')})
    cat=[]
    for c,levels in CATEGORICAL.items():
        for g in [0,1]:
            s=d.loc[d[Y].eq(g),c]; den=int(s.notna().sum())
            for level in levels:
                n=int(s.eq(level).sum())
                cat.append({'variable':c,'IMC':g,'level':level,'n':n,'denominator':den,
                            'percent':100*n/den if den else np.nan})
    missing=pd.DataFrame([{'variable':c,'n':len(d),'missing':int(d[c].isna().sum())} for c in d])
    return pd.DataFrame(rows),pd.DataFrame(cat),missing

def rank_stat(x,labels):
    r=stats.rankdata(np.asarray(x,float),method='average')
    labels=np.asarray(labels,float)
    if labels.ndim==1: labels=labels[None,:]
    n1=labels.sum(axis=1); n0=len(r)-n1
    with np.errstate(divide='ignore',invalid='ignore'):
        out=np.abs((labels@r)/n1-((1-labels)@r)/n0)
    return out

def pearson_stat(codes,labels,k):
    labels=np.asarray(labels,float)
    if labels.ndim==1: labels=labels[None,:]
    H=np.eye(k)[np.asarray(codes,int)]
    totals=H.sum(axis=0); n=len(codes); n1=labels.sum(axis=1)
    obs1=labels@H; obs0=totals-obs1
    exp1=n1[:,None]*totals/n; exp0=(n-n1[:,None])*totals/n
    with np.errstate(divide='ignore',invalid='ignore'):
        a=np.divide((obs1-exp1)**2,exp1,out=np.zeros_like(exp1),where=exp1>0)
        b=np.divide((obs0-exp0)**2,exp0,out=np.zeros_like(exp0),where=exp0>0)
    out=(a+b).sum(axis=1)
    out[(n1==0)|(n1==n)]=np.nan
    return out

def permutation_tests(d,B=N_PERM):
    y=d[Y].to_numpy(float); rng=rng_for('shared_labels')
    labels=np.vstack([y]+[rng.permutation(y) for _ in range(B)])
    rows=[]; distributions={}
    for family,variables in FAMILIES.items():
        for c in variables:
            ok=d[c].notna().to_numpy(); x=d.loc[ok,c]; yy=labels[:,ok]
            if len(x)<3 or x.nunique()<2:
                raise ValueError('Non-testable prespecified variable: '+c)
            if c in CATEGORICAL:
                levels=[v for v in CATEGORICAL[c] if x.eq(v).any()]
                codes=pd.Categorical(x,categories=levels).codes
                t=pearson_stat(codes,yy,len(levels))
            else: t=rank_stat(x.to_numpy(float),yy)
            if not np.isfinite(t).all():
                raise ValueError('Degenerate shared permutation; no silent deletion: '+c)
            tol=1e-10*max(1,abs(t[0]))
            p=(1+np.count_nonzero(t[1:]>=t[0]-tol))/(B+1)
            rows.append({'family':family,'variable':c,'n':len(x),'events':int(y[ok].sum()),
                         'statistic':float(t[0]),'p_permutation':float(p)})
            # Round only for numerical tie recognition, not for reported estimates.
            distributions[c]=stats.rankdata(-np.round(t,10),method='max')/(B+1)
    out=pd.DataFrame(rows); out['q_BH']=np.nan; out['p_WY_minP']=np.nan
    for family,variables in FAMILIES.items():
        ix=out.index[out.family.eq(family)]
        out.loc[ix,'q_BH']=bh(out.loc[ix,'p_permutation'])
        if RUN_SUPPLEMENTARY_WY:
            minp=np.min(np.vstack([distributions[c] for c in variables]),axis=0)
            for i in ix:
                p=out.loc[i,'p_permutation']
                adjusted=(1+np.count_nonzero(minp[1:]<=p+1e-12))/(B+1)
                out.loc[i,'p_WY_minP']=max(p,adjusted)
    return out

def robust_effects(d,B=N_BOOT):
    rows=[]
    for c in CLINICAL_NUM+LAB:
        a=d.loc[d[Y].eq(0),c].dropna().to_numpy(float)
        b=d.loc[d[Y].eq(1),c].dropna().to_numpy(float)
        if min(len(a),len(b))<2: raise ValueError('Insufficient effect-size data: '+c)
        diff=b[:,None]-a[None,:]
        r=float(np.sign(diff).mean()); hl=float(np.median(diff))
        rng=rng_for('bootstrap_'+c); rb=[]; hb=[]
        for start in range(0,B,250):
            size=min(250,B-start)
            aa=rng.choice(a,size=(size,len(a)),replace=True)
            bb=rng.choice(b,size=(size,len(b)),replace=True)
            delta=bb[:,:,None]-aa[:,None,:]
            rb.extend(np.sign(delta).mean(axis=(1,2)))
            hb.extend(np.median(delta,axis=(1,2)))
        rl,rh=np.quantile(rb,[.025,.975]); hl_l,hl_h=np.quantile(hb,[.025,.975])
        rows.append({'variable':c,'n0':len(a),'n1':len(b),'rank_biserial':r,
                     'r_CI_low':rl,'r_CI_high':rh,'HL_shift':hl,'HL_CI_low':hl_l,'HL_CI_high':hl_h,
                     'median_difference':float(np.median(b)-np.median(a))})
    return pd.DataFrame(rows)

def categorical_effects(d):
    rows=[]
    for c,declared in CATEGORICAL.items():
        z=d[[Y,c]].dropna(); levels=[v for v in declared if z[c].eq(v).any()]
        table=pd.crosstab(z[Y],z[c]).reindex(index=[0,1],columns=levels,fill_value=0).to_numpy()
        row={'variable':c,'n':len(z),'levels':len(levels),'V_Cramer':np.nan,
             'OR_conditional':np.nan,'OR_CI_low':np.nan,'OR_CI_high':np.nan,'p_Fisher':np.nan}
        if len(levels)<2 or (table.sum(axis=1)==0).any():
            row['status']='not_testable'; rows.append(row); continue
        chi=stats.chi2_contingency(table,correction=False).statistic
        row['V_Cramer']=float(np.sqrt(chi/len(z)))
        if len(levels)==2:
            # [[IMC second,IMC first],[no-IMC second,no-IMC first]]
            t=table[::-1,::-1]
            fit=odds_ratio(t,kind='conditional'); ci=fit.confidence_interval()
            row.update(OR_conditional=float(fit.statistic),OR_CI_low=float(ci.low),
                       OR_CI_high=float(ci.high),p_Fisher=float(stats.fisher_exact(t).pvalue),
                       comparison=levels[1]+' vs '+levels[0])
        row['status']='ok'; rows.append(row)
    return pd.DataFrame(rows)

def transform_predictor(x,name):
    x=np.asarray(x,float)
    if name in LOG2:
        if (x<=0).any(): raise ValueError('Nonpositive log2 predictor: '+name)
        return np.log2(x),'per doubling',np.nan
    iqr=float(np.quantile(x,.75)-np.quantile(x,.25))
    scale=iqr if iqr>0 else 1.
    return x/scale,('per IQR' if iqr>0 else 'per source unit; IQR=0'),scale

def ml_fit(y,X):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if np.linalg.matrix_rank(X)<X.shape[1]: raise ValueError('rank_deficient')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit=sm.GLM(y,X,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10)
    if not fit.converged or not np.isfinite(fit.params).all(): raise ValueError('ML_not_converged')
    if any('separation' in str(w.message).lower() for w in caught): raise ValueError('ML_separation')
    return fit

def logistic_analyses(d,cat_effect):
    rows=[]
    for family,variables in FAMILIES.items():
        for c in variables:
            z=d[[Y,c]].dropna(); y=z[Y].to_numpy(float)
            row={'family':family,'variable':c,'n':len(z),'events':int(y.sum()),'p':np.nan,
                 'OR':np.nan,'CI_low':np.nan,'CI_high':np.nan,'status':'not_fitted'}
            try:
                if c in CATEGORICAL:
                    levels=[v for v in CATEGORICAL[c] if z[c].eq(v).any()]
                    if len(levels)==2:
                        r=cat_effect.loc[cat_effect.variable.eq(c)].iloc[0]
                        row.update(OR=r.OR_conditional,CI_low=r.OR_CI_low,CI_high=r.OR_CI_high,
                                   p=r.p_Fisher,method='conditional_exact_OR_Fisher',comparison=r.comparison)
                    else:
                        x=pd.get_dummies(pd.Categorical(z[c],categories=levels),drop_first=True,dtype=float)
                        X=np.column_stack([np.ones(len(z)),x.to_numpy()])
                        fit=ml_fit(y,X); null=ml_fit(y,np.ones((len(y),1)))
                        row.update(p=float(stats.chi2.sf(max(0,2*(fit.llf-null.llf)),len(levels)-1)),
                                   method='ML_global_LR',comparison='omnibus')
                else:
                    x,unit,scale=transform_predictor(z[c],c)
                    fit=ml_fit(y,np.column_stack([np.ones(len(y)),x])); ci=fit.conf_int()[1]
                    row.update(OR=float(np.exp(fit.params[1])),CI_low=float(np.exp(ci[0])),
                               CI_high=float(np.exp(ci[1])),p=float(fit.pvalues[1]),
                               method='ML_Wald',comparison=unit,scale=scale)
                row['status']='ok'
            except ValueError as e: row['status']=str(e)
            rows.append(row)
    out=pd.DataFrame(rows); out['q_BH']=np.nan
    for family in FAMILIES:
        ix=out.index[out.family.eq(family)]; out.loc[ix,'q_BH']=bh(out.loc[ix,'p'])
    adjusted=[]
    for label,covs in [('A_age_PICU',['wiek_lata','oitd']),('B_age_WBCmin',['wiek_lata','wbc_min'])]:
        z=d[[Y]+covs].dropna(); cols=[np.ones(len(z))]; units={}
        for c in covs:
            if c=='oitd': cols.append(z[c].map({'nie':0.,'tak':1.}).to_numpy()); units[c]='yes vs no'
            else:
                x,unit,scale=transform_predictor(z[c],c); cols.append(x); units[c]=unit+f' (scale={scale:g})'
        try:
            fit=ml_fit(z[Y],np.column_stack(cols))
            for j,c in enumerate(covs,1):
                ci=fit.conf_int()[j]
                adjusted.append({'model':label,'variable':c,'n':len(z),'events':int(z[Y].sum()),
                    'unit':units[c],'OR':float(np.exp(fit.params[j])),
                    'CI_low':float(np.exp(ci[0])),'CI_high':float(np.exp(ci[1])),
                    'p_Wald':float(fit.pvalues[j]),'status':'ok_descriptive_only'})
        except ValueError as e: adjusted.append({'model':label,'status':str(e)})
    return out,pd.DataFrame(adjusted)

def firth_value_score(beta,X,y):
    eta=X@beta; p=special.expit(eta); w=p*(1-p)
    I=X.T@(w[:,None]*X); sign,ld=np.linalg.slogdet(I)
    if sign<=0: return -np.inf,np.full(len(beta),np.nan)
    try: inv=np.linalg.inv(I)
    except np.linalg.LinAlgError: return -np.inf,np.full(len(beta),np.nan)
    h=w*np.einsum('ij,jk,ik->i',X,inv,X)
    ll=float(np.sum(y*eta-np.logaddexp(0,eta))+.5*ld)
    score=X.T@(y-p+h*(.5-p))
    return ll,score

def firth_opt(X,y,fixed=None,start=None):
    k=X.shape[1]; fixed={} if fixed is None else dict(fixed)
    free=[j for j in range(k) if j not in fixed]
    base=np.zeros(k) if start is None else np.array(start,float).copy()
    for j,v in fixed.items(): base[j]=v
    def objective(v):
        b=base.copy(); b[free]=v; ll,s=firth_value_score(b,X,y)
        return -ll,-s[free]
    res=optimize.minimize(objective,base[free],jac=True,method='BFGS',
                          options={'gtol':1e-8,'maxiter':1000})
    b=base.copy(); b[free]=res.x; ll,s=firth_value_score(b,X,y)
    if not np.isfinite(ll) or (len(free) and np.max(np.abs(s[free]))>2e-5):
        raise ValueError('Firth optimizer failed score convergence.')
    return b,ll

def firth_profile(X,y,j=1):
    b,ll=firth_opt(X,y); cutoff=stats.chi2.ppf(.95,1)
    def target(v):
        _,lp=firth_opt(X,y,{j:float(v)},b)
        return 2*(ll-lp)-cutoff
    def boundary(direction):
        step=.5
        for _ in range(30):
            endpoint=b[j]+direction*step
            if target(endpoint)>0:
                return optimize.brentq(target,min(b[j],endpoint),max(b[j],endpoint),xtol=1e-8)
            step*=1.5
        raise ValueError('Firth profile CI could not be bracketed.')
    low=boundary(-1); high=boundary(1)
    _,lnull=firth_opt(X,y,{j:0.},b)
    p=float(stats.chi2.sf(max(0,2*(ll-lnull)),1))
    return {'OR':float(np.exp(b[j])),'CI_low':float(np.exp(low)),
            'CI_high':float(np.exp(high)),'p_profile':p,'status':'ok'}

def firth_platelet_models(d):
    z=d[[Y,'plt_min']].dropna(); scale=float(z.plt_min.quantile(.75)-z.plt_min.quantile(.25))
    if scale<=0: raise ValueError('Platelet IQR must be positive.')
    rows=[]
    for label,covs in [('unadjusted',[]),('age_adjusted',['wiek_lata'])]:
        z=d[[Y,'plt_min']+covs].dropna(); x=(z.plt_min-z.plt_min.mean())/scale
        X=np.column_stack([np.ones(len(z)),x])
        if covs:
            age=z.wiek_lata.to_numpy(); sd=age.std(ddof=1)
            if sd<=0: raise ValueError('Age is constant.')
            X=np.column_stack([X,(age-age.mean())/sd])
        row={'model':label,'n':len(z),'events':int(z[Y].sum()),'platelet_IQR':scale}
        try: row.update(firth_profile(X,z[Y].to_numpy(float)))
        except ValueError as e: row['status']=str(e)
        rows.append(row)
    return pd.DataFrame(rows)

def full_cohort_outputs(full):
    numeric=[]
    for c in ['wiek_lata']+LAB+[v for v in OPTIONAL_NUM if v in full]:
        x=full[c].dropna()
        numeric.append({'variable':c,'n':len(x),'median':x.median(),'Q1':x.quantile(.25),'Q3':x.quantile(.75)})
    cat=[]
    for c,levels in CATEGORICAL.items():
        s=full[c].dropna()
        if c=='zomr': s=s[s.ne('nie badano')]; levels=['nie','tak']
        for level in levels:
            n=int(s.eq(level).sum()); den=len(s)
            cat.append({'variable':c,'level':level,'n':n,'denominator':den,'percent':100*n/den if den else np.nan})
    year=full['onset_year'].astype(int)
    annual=pd.DataFrame({'year':year,'status':np.where(full.death_within_24h.eq(1),
        'Death within 24 h',np.where(full[Y].eq(1),'With IMC','Without IMC'))})
    annual=pd.crosstab(annual.year,annual.status).reindex(
        index=range(int(year.min()),int(year.max())+1),
        columns=['With IMC','Without IMC','Death within 24 h'],fill_value=0).fillna(0).astype(int)
    return pd.DataFrame(numeric),pd.DataFrame(cat),annual

def annual_figure(annual):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    bottom=np.zeros(len(annual))
    for c in annual.columns:
        ax.bar(annual.index,annual[c],bottom=bottom,label=c); bottom+=annual[c].to_numpy()
    ax.set_xlabel('Year of disease onset'); ax.set_ylabel('Number of cases')
    ax.set_title('Invasive meningococcal disease: annual case counts')
    ax.set_xticks(annual.index); ax.legend(); fig.tight_layout()
    return fig

assert np.isclose(np.sign(np.array([3,4])[:,None]-np.array([1,2])[None,:]).mean(),1)
assert np.allclose(bh([.01,.04,.03]),[.03,.04,.04])
assert np.isclose(rank_stat([1,2,3,4],[0,0,1,1])[0],2)
assert np.isclose(pearson_stat([0,0,1,1],[0,0,1,1],2)[0],4)
assert odds_ratio([[4,1],[2,5]],kind='conditional').statistic>1
Xtest=np.column_stack([np.ones(8),np.arange(8)/4])
ytest=np.array([0,1,0,0,1,0,1,1],float)
btest,lltest=firth_opt(Xtest,ytest)
assert np.max(np.abs(firth_value_score(btest,Xtest,ytest)[1]))<2e-5
print('Software checks passed. No study data loaded.')

RESULTS={}; RUN_METADATA={'status':'not_run_no_data'}; ANNUAL_FIGURE=None
if RUN_PRIVATE_ANALYSIS:
    if not CONFIRM_AUTHORIZED_ENVIRONMENT or not CONFIRM_ANALYSIS_COHORT_ELIGIBILITY:
        raise ValueError('Confirm authorized environment and analysis-cohort eligibility first.')
    raw,input_hash=read_private(DATA_PATH,ALLOW_COLAB_UPLOAD)
    d,audit=prepare(raw); del raw
    numeric,cat,missing=descriptive(d)
    RESULTS.update(audit=audit,descriptive_numeric=numeric,descriptive_categorical=cat,missingness=missing)
    RESULTS['permutation_tests']=permutation_tests(d)
    RESULTS['source_discrepancy_flags']=pd.DataFrame([{
        'variable':'przyjecie_do_atb_dni',
        'issue':'Source HTML and manuscript Table 2 disagree; reconcile before publication.',
        'recomputed_p':float(RESULTS['permutation_tests'].loc[
            RESULTS['permutation_tests'].variable.eq('przyjecie_do_atb_dni'),'p_permutation'].iloc[0])}])
    RESULTS['robust_effects']=robust_effects(d)
    RESULTS['categorical_effects']=categorical_effects(d)
    RESULTS['univariable_associations'],RESULTS['adjusted_descriptive_models']=logistic_analyses(d,RESULTS['categorical_effects'])
    if RUN_FIRTH: RESULTS['Firth_platelet_models']=firth_platelet_models(d)
    RESULTS['table2_numeric']=RESULTS['permutation_tests'].merge(RESULTS['robust_effects'],on='variable',how='inner')
    RESULTS['table2_categorical']=RESULTS['permutation_tests'].merge(RESULTS['categorical_effects'],on='variable',how='inner')
    if FULL_COHORT_PATH is not None:
        raw_full,_=read_private(FULL_COHORT_PATH); full,full_audit=prepare(raw_full,full_cohort=True); del raw_full
        if len(full)<=len(d): raise ValueError('Full cohort must include additional unassessable patients.')
        # Aggregate eligibility consistency; private record linkage remains the authors responsibility.
        if int(full[Y].notna().sum())!=len(d): raise ValueError('Full-cohort assessable count differs from analysis cohort.')
        t1n,t1c,annual=full_cohort_outputs(full)
        RESULTS.update(full_cohort_numeric=t1n,full_cohort_categorical=t1c,annual_counts=annual.reset_index())
        ANNUAL_FIGURE=annual_figure(annual)
        import matplotlib.pyplot as plt
        plt.close(ANNUAL_FIGURE)
        del full
    RUN_METADATA={'status':'private_analysis_completed','utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256_PRIVATE':input_hash,'python':platform.python_version(),
        'packages':{p:version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib','nbformat']},
        'seed':SEED,'permutations':N_PERM,'bootstrap':N_BOOT,
        'families':FAMILIES,'missingness':'available cases per variable/model; no imputation',
        'historical_random_stream_recovered':False,'full_cohort_supplied':FULL_COHORT_PATH is not None,
        'timing_verified_from_measurement_timestamps':False}
    if EXPORT_PRIVATE_AGGREGATES:
        PRIVATE_OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
        for name,table in RESULTS.items(): table.to_csv(PRIVATE_OUTPUT_DIR/(name+'.csv'),index=False)
        (PRIVATE_OUTPUT_DIR/'run_metadata_PRIVATE.json').write_text(json.dumps(RUN_METADATA,indent=2),encoding='utf-8')
        if ANNUAL_FIGURE is not None: ANNUAL_FIGURE.savefig(PRIVATE_OUTPUT_DIR/'annual_counts_PRIVATE.png',dpi=300)
    print('Private analysis complete. Aggregate results are in RESULTS; no participant rows exported.')
else:
    print('Public mode: no upload, no study data, no research results computed.')

CLEAN_SOURCE_PATH=None
CLEAN_TARGET_PATH='IMD_IMC_analysis_clean.ipynb'
if CLEAN_SOURCE_PATH is not None:
    import nbformat
    clean=nbformat.read(CLEAN_SOURCE_PATH,as_version=4)
    clean.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},
                    'language_info':{'name':'python'},'colab':{'private_outputs':True,'provenance':[]}}
    for cell in clean.cells:
        cell.metadata={}; cell.pop('attachments',None)
        if cell.cell_type=='code': cell.outputs=[]; cell.execution_count=None
    nbformat.validate(clean); nbformat.write(clean,CLEAN_TARGET_PATH)
    print('Outputs stripped. Review all source text and public defaults before release.')
