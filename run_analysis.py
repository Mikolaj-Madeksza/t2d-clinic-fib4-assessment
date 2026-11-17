"""Run study analyses on artificial examples only; outputs are not study results."""
from derive import (OUT, ID, PID, np, pd, json, derive, select, score_summary,
                    wilson, proportion, continuous, table_characteristics,
                    table_scores, write_json)
from mock_data import make_tables

def main():
    OUT.mkdir(exist_ok=True)
    tables = make_tables()
    e = derive(tables)
    index, primary, flow = select(e)
    pd.DataFrame(flow).to_csv(OUT/'cohort_flow.csv',index=False)
    summaries = {'primary':score_summary(primary), 'by_age':{name:score_summary(g) for name,g in primary.groupby(np.where(primary.age.gt(65),'>65','35–65'))}}
    summaries['availability'] = {'n':len(index),'same_day':int(index.has_panel.sum()),'same_day_ci':wilson(index.has_panel.sum(),len(index)),
                                 'ast':int(index.ast_available.sum()),'alt':int(index.alt_available.sum()),'platelets':int(index.platelets_available.sum()),
                                 'all_components':int(index.all_components.sum())}
    panel = index[index.has_panel]
    summaries['index'] = {'age':index.age_admission.quantile([.25,.5,.75]).tolist(), 'infection':panel.infection.value_counts().to_dict(),
                          'infection_reasons':panel.fib4_infection_reason.value_counts().to_dict(),
                          'cirrhosis':index.cirrhosis_status_at_encounter.value_counts().to_dict(),'pbc':int(index.pbc_documented.sum()),
                          'bmi_source':index.bmi_source.value_counts().to_dict(),
                          'panel_days':(panel.day-panel.admission_date).dt.days.value_counts().sort_index().to_dict()}
    summaries['primary_cirrhosis'] = primary.cirrhosis_status_at_encounter.value_counts().to_dict()
    summaries['primary_crp'] = {'available':int(primary.crp.notna().sum()),'low':int(primary.crp_low.sum()),
                                'censored_low':int((primary.crp_low&primary.crp_operator.notna()).sum())}
    scenarios = [('Primary',index,primary)]
    scenarios += [('CRP ≤1 mg/dL',index,primary[primary.crp_low]),
                  ('Excluding uncertain infection status',index,primary[primary.infection.eq('false')]),
                  ('No suspected/indeterminate cirrhosis',index,primary[primary.cirrhosis_status_at_encounter.eq('no_documented_cirrhosis_after_review')])]
    for name, kwargs in [('All admission modes',{'all_modes':True}),('Youngest plausible age',{'month':12,'day':31}),
                         ('Oldest plausible age',{'month':1,'day':1}),('Higher cutoff from age 65',{'cutoff65':True})]:
        ix,d,_ = select(e,**kwargs);scenarios.append((name,ix,d))
    sensitivity, changes = [], []
    summaries['sensitivities']={}
    for name, ix, d in scenarios:
        summaries['sensitivities'][name]=score_summary(d)
        sensitivity.append([name,str(len(ix)),str(int(ix.has_panel.sum())),str(len(d)),proportion(d.positive.sum(),len(d)),
                            '–'.join(f'{x:.1f}' for x in wilson(d.positive.sum(),len(d))),proportion(d.high.sum(),len(d))])
        match=primary[[PID,ID,'positive','category']].merge(d[[PID,ID,'positive','category']],on=PID,suffixes=('_p','_s'))
        same=match[ID+'_p'].eq(match[ID+'_s'])
        changes.append([name,str(len(match)),str((~same).sum()),str(len(d)-len(match)),str(len(primary)-len(match)),
                        str((same&match.positive_p.ne(match.positive_s)).sum()),str((same&match.category_p.ne(match.category_s)).sum())])
    # Missing components are described separately from a lack of temporal overlap.
    missing=[]
    for key,g in index.groupby(['ast_available','alt_available','platelets_available'],sort=False):
        missing.append({'AST':bool(key[0]),'ALT':bool(key[1]),'platelets':bool(key[2]),'patients':len(g),'same_day':int(g.has_panel.sum())})
    pd.DataFrame(missing).to_csv(OUT/'component_patterns.csv',index=False)
    av=[]
    for (clinic,year),g in index.groupby(['clinic','year']):
        av.append([clinic,str(year),str(len(g)),proportion(g.ast_available.sum(),len(g)),proportion(g.alt_available.sum(),len(g)),
                   proportion(g.platelets_available.sum(),len(g)),proportion(g.has_panel.sum(),len(g))])
    groups=[('Overall',primary)]+list(primary.groupby('clinic'))+list(primary.groupby('period'))
    strat=[]
    for name,g in groups:
        strat.append([name,str(len(g)),continuous(g.fib4,2),proportion(g.positive.sum(),len(g)),
                      '–'.join(f'{x:.1f}' for x in wilson(g.positive.sum(),len(g))),proportion(g.high.sum(),len(g))])
    review=[]
    for status in ['true','false','unknown']:
        review.append(['Acute infection at panel: '+status,proportion(panel.infection.eq(status).sum(),len(panel))])
    for reason,n in panel[panel.infection.eq('unknown')].fib4_infection_reason.value_counts().items():
        label={'infection_timing_or_nature_unresolved':'Unresolved nature/timing','antimicrobial_indication_or_timing_unresolved':'Unresolved antimicrobial context','insufficient_source_context':'Insufficient source context'}[reason]
        review.append(['Unknown infection: '+label,proportion(n,len(panel))])
    for status in ['no_documented_cirrhosis_after_review','established','suspected','indeterminate']:
        review.append(['Cirrhosis: '+status.replace('_',' '),proportion(index.cirrhosis_status_at_encounter.eq(status).sum(),len(index))])
    review.append(['PBC documented',proportion(index.pbc_documented.sum(),len(index))])
    for label,mask in [('Primary cohort: infection unknown',primary.infection.eq('unknown')),
                       ('Primary cohort: suspected cirrhosis',primary.cirrhosis_status_at_encounter.eq('suspected')),
                       ('Primary cohort: indeterminate cirrhosis',primary.cirrhosis_status_at_encounter.eq('indeterminate'))]:
        review.append([label,proportion(mask.sum(),len(primary))])
    publication={'table1':table_characteristics(index),'table2':table_scores(primary),
                 's2':{'columns':['Clinic','Year','Patients','AST','ALT','Platelets','Same-day panel'],'rows':av},
                 's3':{'columns':['Classification','n/N (%)'],'rows':review},
                 's4':{'columns':['Group','Patients','FIB-4, median [IQR]','Threshold-positive','95% CI, %','FIB-4 >2.67'],'rows':strat},
                 's5':{'columns':['Analysis','Index patients','Same-day panels','FIB-4 patients','Threshold-positive','95% CI, %','FIB-4 >2.67'],'rows':sensitivity},
                 's6':{'columns':['Analysis','Common patients','Different index','Added patients','Removed patients','Threshold changes¹','Category changes¹'],'rows':changes}}
    summaries['flow']=flow
    summaries['threshold_difference_pp']=100*(primary.fixed_positive.sum()-primary.positive.sum())/len(primary)
    write_json(OUT/'summary.json',summaries);write_json(OUT/'tables.json',publication)
    for name,t in publication.items():pd.DataFrame(t['rows'],columns=t['columns']).to_csv(OUT/(name+'.csv'),index=False)
    print(json.dumps({k:summaries[k] for k in ['availability','primary','primary_crp','primary_cirrhosis','flow']},indent=2,default=str))

if __name__ == '__main__':
    main()
