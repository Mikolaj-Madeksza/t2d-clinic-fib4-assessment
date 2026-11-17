"""Study derivations adapted without clinical file loading or patient-level exports."""
from pathlib import Path
import json
from statistics import NormalDist
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / 'outputs'
ID, PID = 'encounter_study_id', 'patient_study_id'
CORE = {'ast_u_l': ('ast', 'U/L'), 'alt_u_l': ('alt', 'U/L'),
        'platelets_10e3_ul': ('platelets', '10^3/uL')}
UNITS = {a: u for a, (_, u) in CORE.items()} | {'hba1c_percent': '%', 'crp_mg_dl': 'mg/dL'}


def age_at(date, birth_year, month=7, day=1):
    return date.dt.year - birth_year - ((date.dt.month * 100 + date.dt.day) < month * 100 + day).astype(int)


def wilson(k, n):
    if not n:
        return [None, None]
    z = NormalDist().inv_cdf(.975)
    p = k / n
    a = (p + z*z/(2*n)) / (1 + z*z/n)
    b = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n)
    return [100*(a-b), 100*(a+b)]


def proportion(k, n, ci=False):
    if not n:
        return '—'
    text = f'{int(k)}/{int(n)} ({100*k/n:.1f}%)'
    if ci:
        lo, hi = wilson(k, n)
        text += f'; {lo:.1f}–{hi:.1f}'
    return text


def continuous(values, digits=1):
    v = pd.Series(values).dropna()
    if not len(v):
        return '—; n=0'
    q1, med, q3 = v.quantile([.25, .5, .75])
    return f'{med:.{digits}f} [{q1:.{digits}f}–{q3:.{digits}f}]; n={len(v)}'


def derive(t):
    e = t['encounters'].copy()
    for col in ['admission_date', 'discharge_date']:
        e[col] = pd.to_datetime(e[col])
    e['birth_year'] = e.admission_date.dt.year - e.age_at_admission_calendar_years
    assert e.groupby(PID).birth_year.nunique().eq(1).all()
    assert e.groupby(PID).sex.nunique().eq(1).all()
    h = t['clinical_history']
    pregnant = h[h.concept.eq('pregnancy') & h.assertion.eq('present') & ~h.temporality.eq('historical')]
    e['pregnancy'] = e[ID].isin(pregnant[ID])
    e['year'] = e.admission_date.dt.year
    e['period'] = np.where(e.year.le(2022), '2017–2022', '2023–2026')

    q = t['labs'][t['labs'].analyte.isin(UNITS)].copy()
    q = q.merge(e[[ID, PID, 'admission_date', 'discharge_date']], on=[ID, PID], how='left', validate='many_to_one')
    q['day'] = pd.to_datetime(q.collected_at).dt.normalize()
    q['usable'] = (q.encounter_link_status.eq('linked') & q.day.between(q.admission_date, q.discharge_date)
                   & q.result_status.eq('parsed') & np.isfinite(q.value) & q.unit.eq(q.analyte.map(UNITS)))
    q['exact'] = q.usable & q.comparison_operator.isna() & q.value.gt(0)
    records = []
    for a, g in q.groupby('analyte'):
        records.append({'analyte': a, 'all_records': len(g), 'missing_collection': int(g.day.isna().sum()),
                        'nonexact_records': int(g.comparison_operator.notna().sum()),
                        'valid_exact_in_stay': int(g.exact.sum())})
    pd.DataFrame(records).to_csv(OUT / 'laboratory_checks.csv', index=False)
    good = q[q.exact].sort_values(['collected_at', 'lab_result_id'])
    daily = good[good.analyte.isin(CORE)].drop_duplicates([ID, 'day', 'analyte'])
    complete = daily.pivot(index=[ID, 'day'], columns='analyte', values='value').dropna(subset=list(CORE)).reset_index()
    first = complete.sort_values([ID, 'day']).drop_duplicates(ID)
    first = first.rename(columns={a: name for a, (name, _) in CORE.items()})
    e = e.merge(first, on=ID, how='left', validate='one_to_one')
    e['has_panel'] = e.day.notna()
    for a, (name, _) in CORE.items():
        earliest = good[good.analyte.eq(a)].drop_duplicates(ID).set_index(ID).value
        e[f'{name}_available'] = e[ID].isin(earliest.index)
        e[f'{name}_descriptive'] = e[name].fillna(e[ID].map(earliest))
    e['all_components'] = e[['ast_available', 'alt_available', 'platelets_available']].all(axis=1)
    e['enzyme_ok'] = e.has_panel & e.ast.le(200) & e.alt.le(200)
    e['infection'] = np.select([e.acute_infection_at_fib4_panel.eq(True), e.acute_infection_at_fib4_panel.eq(False)], ['true', 'false'], default='unknown')
    e.loc[~e.has_panel, 'infection'] = 'not_assessable'
    assert e.loc[e.has_panel, 'fib4_infection_review_status'].isin(['completed', 'completed_unknown']).all()
    assert e.loc[~e.has_panel, 'fib4_infection_review_status'].eq('not_assessable_no_panel').all()
    assert e.loc[e.infection.eq('unknown'), 'fib4_infection_review_status'].eq('completed_unknown').all()
    assert e.loc[e.infection.isin(['true', 'false']), 'fib4_infection_review_status'].eq('completed').all()

    hb = good[good.analyte.eq('hba1c_percent')].drop_duplicates(ID).set_index(ID).value
    e['hba1c'] = e[ID].map(hb)
    crp = q[q.usable & q.analyte.eq('crp_mg_dl') & q.value.ge(0)
            & (q.comparison_operator.isna() | q.comparison_operator.isin(['<', '<=']))].sort_values(['collected_at', 'lab_result_id'])
    cs = crp.drop_duplicates([ID, 'day'])[[ID, 'day', 'value', 'comparison_operator']]
    e = e.merge(cs.rename(columns={'value': 'crp', 'comparison_operator': 'crp_operator'}), on=[ID, 'day'], how='left', validate='one_to_one')
    ce = crp.drop_duplicates(ID).set_index(ID)
    e['crp_descriptive'] = e.crp.where(e.has_panel, e[ID].map(ce.value))
    e['crp_descriptive_operator'] = e.crp_operator.where(e.has_panel, e[ID].map(ce.comparison_operator))
    e['crp_exact'] = e.crp_descriptive.where(e.crp_descriptive_operator.isna())
    e['crp_censored'] = e.crp_descriptive.notna() & e.crp_descriptive_operator.notna()
    e['crp_low'] = e.crp.notna() & e.crp.le(1)

    m = t['measurements']
    bm = m[m.measurement.eq('body_mass_index') & m.value_quality_status.eq('plausible')
           & m.unit.eq('kg/m²') & m.measurement_within_encounter_dates.eq(True)
           & m.value_numeric.gt(0) & np.isfinite(m.value_numeric) & m[ID].notna()]
    means = bm.groupby(ID).value_numeric.mean()
    e['bmi'] = e[ID].map(means).fillna(e.calculated_bmi_kg_m2).fillna(e.documented_bmi_kg_m2)
    assert e.bmi.dropna().gt(0).all()
    e['bmi_source'] = np.select([e[ID].isin(means.index), e.calculated_bmi_kg_m2.notna(), e.documented_bmi_kg_m2.notna()],
                              ['mean_in_stay', 'calculated_admission', 'documented_admission'], default='missing')
    return e


def select(e, month=7, day=1, all_modes=False, cutoff65=False):
    e = e.copy()
    e['age_admission'] = age_at(e.admission_date, e.birth_year, month, day)
    e['age'] = age_at(e.day, e.birth_year, month, day)
    flow = []
    def step(label, mask):
        nonlocal e
        before = len(e)
        e = e[mask.reindex(e.index)].copy()
        flow.append({'step': label, 'admissions': len(e), 'patients': e[PID].nunique(), 'excluded_admissions': before-len(e)})
    step('Reviewed T2D input cohort', pd.Series(True, index=e.index))
    step('Discharged adults in study dates', e.hospitalization_status.eq('discharged') & e.age_admission.ge(18) & e.admission_date.between('2017-01-01', '2026-09-03'))
    step('No current pregnancy', ~e.pregnancy)
    step('Age at least 35 years', e.age_admission.ge(35))
    if not all_modes:
        step('Planned admission', e.admission_mode_group.eq('planned'))
    step('No established cirrhosis', ~e.cirrhosis_status_at_encounter.eq('established'))
    index = e.sort_values(['admission_date', ID]).drop_duplicates(PID).copy()
    flow.append({'step': 'Index patients', 'admissions': len(index), 'patients': len(index), 'excluded_admissions': len(e)-len(index)})
    index['lower_cutoff'] = np.where(index.age.ge(65) if cutoff65 else index.age.gt(65), 2., 1.3)
    index['fib4'] = index.age * index.ast / (index.platelets * np.sqrt(index.alt))
    index['positive'] = index.fib4.ge(index.lower_cutoff)
    index['fixed_positive'] = index.fib4.ge(1.3)
    index['high'] = index.fib4.gt(2.67)
    index['category'] = np.select([index.high, index.positive], ['Higher', 'Intermediate'], default='Lower')
    index.loc[~index.has_panel, 'category'] = None
    e = index.copy()
    step('Valid same-day panel', e.has_panel)
    step('AST and ALT at most 200 U/L', e.enzyme_ok)
    step('Primary FIB-4 cohort', ~e.infection.eq('true'))
    assert e[PID].is_unique and e.infection.isin(['false', 'unknown']).all()
    return index, e, flow


def score_summary(d):
    n = len(d)
    return {'n': n, 'positive': int(d.positive.sum()), 'positive_pct': 100*d.positive.mean() if n else None,
            'positive_ci': wilson(int(d.positive.sum()), n), 'high': int(d.high.sum()),
            'high_ci': wilson(int(d.high.sum()), n), 'fixed': int(d.fixed_positive.sum()),
            'lower': int(d.category.eq('Lower').sum()), 'intermediate': int(d.category.eq('Intermediate').sum()),
            'reclassified': int((d.fixed_positive & ~d.positive).sum()),
            'fib4_quantiles': d.fib4.quantile([.25,.5,.75]).tolist(), 'infection_unknown': int(d.infection.eq('unknown').sum())}


def table_characteristics(index):
    groups = [('All index patients', index), ('Same-day panel', index[index.has_panel]), ('No same-day panel', index[~index.has_panel])]
    rows = []
    definitions = [('Patients', lambda d: str(len(d))), ('Age, years', lambda d: continuous(d.age_admission)),
                   ('Female', lambda d: proportion(d.sex.eq('female').sum(), len(d))),
                   ('Diabetology', lambda d: proportion(d.clinic.eq('Diabetology').sum(), len(d))),
                   ('Endocrinology', lambda d: proportion(d.clinic.eq('Endocrinology').sum(), len(d))),
                   ('Admission year', lambda d: continuous(d.year, 0)),
                   ('BMI, kg/m²', lambda d: continuous(d.bmi)),
                   ('BMI unavailable', lambda d: proportion(d.bmi.isna().sum(), len(d))),
                   ('HbA1c, %', lambda d: continuous(d.hba1c)),
                   ('HbA1c unavailable', lambda d: proportion(d.hba1c.isna().sum(), len(d))),
                   ('AST, U/L', lambda d: continuous(d.ast_descriptive, 0)),
                   ('ALT, U/L', lambda d: continuous(d.alt_descriptive, 0)),
                   ('Platelets, 10⁹/L', lambda d: continuous(d.platelets_descriptive, 0)),
                   ('CRP, mg/dL (exact results)', lambda d: continuous(d.crp_exact, 2)),
                   ('CRP below reporting limit', lambda d: proportion(d.crp_censored.sum(), len(d))),
                   ('CRP unavailable', lambda d: proportion(d.crp_descriptive.isna().sum(), len(d))),
                   ('AST available', lambda d: proportion(d.ast_available.sum(), len(d))),
                   ('ALT available', lambda d: proportion(d.alt_available.sum(), len(d))),
                   ('Platelets available', lambda d: proportion(d.platelets_available.sum(), len(d))),
                   ('All three components available', lambda d: proportion(d.all_components.sum(), len(d))),
                   ('Complete same-day panel', lambda d: proportion(d.has_panel.sum(), len(d)))]
    for label, fn in definitions:
        rows.append([label] + [fn(d) for _,d in groups])
    return {'columns': ['Characteristic']+[f'{name}\nN={len(d)}' for name,d in groups], 'rows': rows}


def table_scores(primary):
    groups = [('Overall',primary), ('35–65 years',primary[primary.age.le(65)]), ('>65 years',primary[primary.age.gt(65)])]
    rows=[]
    for name, fn in [('Patients',lambda d:str(len(d))), ('FIB-4, median [IQR]',lambda d:continuous(d.fib4,2)),
                     ('Lower category',lambda d:proportion(d.category.eq('Lower').sum(),len(d))),
                     ('Intermediate category',lambda d:proportion(d.category.eq('Intermediate').sum(),len(d))),
                     ('Higher category (>2.67)',lambda d:proportion(d.high.sum(),len(d))),
                     ('Age-adapted threshold-positive',lambda d:proportion(d.positive.sum(),len(d))),
                     ('95% CI, %',lambda d:'–'.join(('—' if v is None else f'{v:.1f}') for v in wilson(d.positive.sum(),len(d)))),
                     ('FIB-4 >2.67: 95% CI, %',lambda d:'–'.join(('—' if v is None else f'{v:.1f}') for v in wilson(d.high.sum(),len(d)))),
                     ('Fixed threshold ≥1.3',lambda d:proportion(d.fixed_positive.sum(),len(d))),
                     ('Reclassified below threshold',lambda d:proportion((d.fixed_positive&~d.positive).sum(),len(d)))]:
        rows.append([name]+[fn(d) for _,d in groups])
    return {'columns':['Measure']+[x for x,_ in groups], 'rows':rows}


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=lambda x: int(x) if isinstance(x,np.integer) else float(x) if isinstance(x,np.floating) else str(x)))
