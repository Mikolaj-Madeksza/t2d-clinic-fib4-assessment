"""Independent checks using artificial records only."""
import math
import numpy as np
import pandas as pd
from mock_data import make_tables
from derive import OUT, ID, PID, derive, select, wilson


def main():
    OUT.mkdir(exist_ok=True)
    tables = make_tables()
    e = derive(tables)
    index, primary, _ = select(e)
    # Independently select first eligible admissions from raw artificial rows.
    pregnant = {r[ID] for r in tables['clinical_history'].to_dict('records')
                if r['assertion'] == 'present' and r['temporality'] != 'historical'}
    expected = {}
    for r in sorted(tables['encounters'].to_dict('records'), key=lambda x: (x['admission_date'], x[ID])):
        age = r['age_at_admission_calendar_years'] - (r['admission_date'].month < 7)
        if (r['hospitalization_status'] == 'discharged' and age >= 35
                and r['admission_mode_group'] == 'planned'
                and r['cirrhosis_status_at_encounter'] != 'established' and r[ID] not in pregnant):
            expected.setdefault(r[PID], r[ID])
    assert dict(zip(index[PID], index[ID])) == expected
    for i in [10, 11, 12, 13]:
        encounter = f'SYNTHETIC-ADMISSION-{i:03}'
        assert encounter in set(index[ID]) and encounter not in set(primary[ID])
    for row in primary.to_dict('records'):
        raw = tables['labs']
        panel = raw[(raw[ID] == row[ID]) & raw.analyte.isin(['ast_u_l', 'alt_u_l', 'platelets_10e3_ul'])]
        values = dict(zip(panel.analyte, panel.value))
        score = row['age'] * values['ast_u_l'] / (values['platelets_10e3_ul'] * math.sqrt(values['alt_u_l']))
        assert math.isclose(score, row['fib4'], rel_tol=1e-12)
        assert row['positive'] == (score >= (2 if row['age'] > 65 else 1.3))
        assert row['high'] == (score > 2.67)
    assert set(primary[primary.infection.eq('unknown')][ID]) == {'SYNTHETIC-ADMISSION-014'}
    assert not index.loc[index[ID].eq('SYNTHETIC-ADMISSION-011'), 'has_panel'].iloc[0]
    assert index.loc[index[ID].eq('SYNTHETIC-ADMISSION-011'), 'all_components'].iloc[0]
    assert e.loc[e[ID].eq('SYNTHETIC-ADMISSION-020'), 'crp_low'].iloc[0]
    # Invalid first AST must not be usable; later valid panels remain selectable.
    for field, value in [('unit', 'wrong'), ('comparison_operator', '<'), ('result_status', 'invalid'), ('value', -1.)]:
        modified = {k: v.copy() for k, v in tables.items()}
        mask = modified['labs'][ID].eq('SYNTHETIC-ADMISSION-020') & modified['labs'].analyte.eq('ast_u_l')
        modified['labs'].loc[mask, field] = value
        modified['encounters'].loc[modified['encounters'][ID].eq('SYNTHETIC-ADMISSION-020'), 'fib4_infection_review_status'] = 'not_assessable_no_panel'
        altered = derive(modified)
        assert not altered.loc[altered[ID].eq('SYNTHETIC-ADMISSION-020'), 'has_panel'].iloc[0]
    for kwargs in [{'all_modes': True}, {'month': 12, 'day': 31}, {'month': 1, 'day': 1}, {'cutoff65': True}]:
        ix, cohort, _ = select(e, **kwargs)
        assert ix[PID].is_unique and set(cohort[ID]) <= set(ix[ID])
        assert cohort.has_panel.all() and cohort.enzyme_ok.all() and not cohort.infection.eq('true').any()
    lo, hi = wilson(5, 10)
    assert math.isclose(lo, 23.6593090513, abs_tol=1e-8) and math.isclose(hi, 76.3406909487, abs_tol=1e-8)
    assert wilson(0, 0) == [None, None]
    assert all(t[ID].str.startswith('SYNTHETIC-').all() and t[PID].str.startswith('SYNTHETIC-').all() for t in tables.values())
    print('Synthetic selection, calculations, missingness, invalid-lab and sensitivity checks passed.')


if __name__ == '__main__':
    main()
