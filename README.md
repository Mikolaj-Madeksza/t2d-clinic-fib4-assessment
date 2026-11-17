# FIB-4 assessment in adults with type 2 diabetes

Code for a descriptive study of age-adapted FIB-4 categories at planned specialist
hospital admissions, with laboratory availability as a supporting outcome.

## Data and scope

**No real patient data or study results are included.** `mock_data.py` creates
four linked tables in memory using wholly invented values. Examples are not
sampled, perturbed or reconstructed from patients, and are not fitted to study
distributions. All identifiers begin with `SYNTHETIC-`. The scripts read no
clinical files. Synthetic outputs demonstrate execution, not the paper's results.

This package follows the companion kidney-testing and dyslipidemia repositories.
The derivations and descriptive analyses are adapted from the performed study
analysis; clinical file loading, source-release verification, patient-level
exports, manuscripts and figures are omitted. Upstream clinical review is
represented by documented input classifications, not recreated by this code.

## Files and use

- `derive.py`: laboratory selection, eligibility, FIB-4 and descriptive summaries.
- `run_analysis.py`: cohort flow, main tables and supporting/sensitivity analyses.
- `mock_data.py`: artificial input tables, including repeated admissions and exclusions.
- `validate.py`: independent selection and score checks, plus edge cases.
- `variable_dictionary.csv`: all input fields and principal derived variables.
- `code_lists.csv`: clinical concepts, units and selection rules.

Use Python 3.14 and the pinned dependencies:

```sh
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_analysis.py
python validate.py
```

Only aggregate artificial outputs are written to `outputs/`; records remain in
memory. Outputs and environments are excluded from Git. The `.gitignore` permits
only the ten reviewed source/documentation files; do not force-add clinical data.

## Selection and measurement

The source cohort was identified by ICD-10 E11 in any diagnosis position, with
type 2 diabetes verified from clinical documentation. That review precedes these
inputs. Hospital tables are linked by patient and admission identifiers.

Select discharged adults admitted from 1 January 2017 through 3 September 2026,
exclude current pregnancy and age below 35, restrict to planned admissions, and
exclude established cirrhosis. Select each patient's first eligible admission
**before** assessing laboratory availability, values or infection. A later
admission never substitutes for an unassessable index admission.

FIB-4 uses the earliest complete same-day AST, ALT and platelet panel during the
index admission. Components must be exact positive values in the required units;
collection time and laboratory identifier resolve duplicate results. Reject
selected panels with AST or ALT above 200 U/L without substituting later panels.
Documented acute infection active at sampling excludes score assessment;
uncertain status remains in the primary analysis. No documented infection is not
proof of clinical absence. Infection and cirrhosis classifications depend on the
available clinical documentation; later diagnoses are not backdated.

Age assumes a 1 July birth date. Admission age determines eligibility; sampling
age enters FIB-4 = age × AST / (platelets × √ALT). Lower thresholds are 1.3 at
35–65 years and 2.0 above 65. Intermediate scores extend from the lower threshold
through 2.67; higher scores exceed 2.67. These are assessment categories, not
confirmed fibrosis or completed referrals. Platelets in 10³/µL equal 10⁹/L.

## Analyses

Availability uses all index patients. FIB-4 categories use eligible index panels.
Summaries report n/N (%), medians and interquartile ranges, and Wilson 95%
confidence intervals. The fixed 1.3 threshold and age-adapted threshold are
compared in the same patients. Supporting outputs describe availability by clinic
and year, missing components, clinical classifications and score strata.

Sensitivity analyses exclude uncertain infection, exclude suspected/indeterminate
cirrhosis, restrict to same-day CRP ≤1 mg/dL, include all admission modes, vary
birth dates to 1 January/31 December, and apply the higher lower threshold from
age 65. CRP below a qualifying reporting limit is eligible without assigning an
exact concentration. Age and admission-mode analyses repeat index selection;
other restrictions retain the original admissions. Changes in cohort membership
are distinguished from classification changes at the same admission.

There are no regression models, hypothesis tests or imputations. Laboratory
unavailability does not prove a test was not performed. Partial coverage in 2026
and changes in clinic organization limit temporal interpretation. Clinical
extraction/classification accuracy is not established by these code checks.
Code and artificial examples use the MIT License.
