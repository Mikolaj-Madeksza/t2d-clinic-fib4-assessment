"""Wholly invented fixtures, unrelated to any patient or study distribution.
No external files are read. Identifiers and dates are generated here.
"""
import numpy as np
import pandas as pd

ID, PID = 'encounter_study_id', 'patient_study_id'


def make_tables():
    tables = {k: [] for k in ['encounters', 'labs', 'clinical_history', 'measurements']}
    for i in range(80):
        person = i % 40
        date = pd.Timestamp(2023 if i < 40 else 2025, 8, 1)
        birth = 2023 - (30 + person)
        keys = {ID: f'SYNTHETIC-ADMISSION-{i:03}', PID: f'SYNTHETIC-PATIENT-{person:03}'}
        missing = i == 10
        split = i == 11
        infection = True if i == 13 else None if i == 14 else False
        tables['encounters'].append({**keys, 'admission_date': date,
            'discharge_date': date + pd.Timedelta(days=4),
            'age_at_admission_calendar_years': date.year - birth,
            'sex': 'female' if person % 2 else 'male',
            'clinic': 'Diabetology' if person % 2 else 'Endocrinology',
            'hospitalization_status': 'ongoing' if i == 6 else 'discharged',
            'admission_mode_group': 'emergency' if i == 7 else 'planned',
            'cirrhosis_status_at_encounter': {8: 'established', 15: 'suspected', 16: 'indeterminate'}.get(i, 'no_documented_cirrhosis_after_review'),
            'acute_infection_at_fib4_panel': infection,
            'fib4_infection_review_status': 'not_assessable_no_panel' if missing or split else 'completed_unknown' if infection is None else 'completed',
            'fib4_infection_reason': 'insufficient_source_context' if infection is None else 'not_applicable',
            'pbc_documented': i == 17,
            'calculated_bmi_kg_m2': np.nan if i % 3 == 0 else 26.,
            'documented_bmi_kg_m2': 28. if i % 5 == 0 else np.nan})
        tables['clinical_history'].append({**keys, 'concept': 'pregnancy',
            'assertion': 'present' if i == 9 else 'explicitly_absent', 'temporality': 'current'})
        tables['measurements'].append({**keys, 'measurement': 'body_mass_index',
            'value_quality_status': 'plausible', 'unit': 'kg/m²',
            'measurement_within_encounter_dates': True, 'value_numeric': 25. if i % 2 else np.nan})
        values = [('ast_u_l', 220. if i == 12 else 18. + person, 'U/L'),
                  ('alt_u_l', 25., 'U/L'), ('platelets_10e3_ul', 130. + 4*person, '10^3/uL'),
                  ('hba1c_percent', 7., '%'), ('crp_mg_dl', .5 if i % 3 else 2., 'mg/dL')]
        for j, (analyte, value, unit) in enumerate(values):
            if missing and j == 2:
                continue
            tables['labs'].append({**keys, 'lab_result_id': f'SYNTHETIC-LAB-{i:03}-{j}',
                'analyte': analyte, 'value': value, 'unit': unit,
                'collected_at': date + pd.Timedelta(days=1 if split and j == 2 else 0, hours=8),
                'encounter_link_status': 'linked', 'result_status': 'parsed',
                'comparison_operator': '<' if j == 4 and i % 4 == 0 else None})
    return {k: pd.DataFrame(v) for k, v in tables.items()}
