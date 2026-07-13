#!/usr/bin/env bash

set -euo pipefail

# ==========================================================================
# User-editable Hyak / Klone paths and allocation
# ==========================================================================
REPO_DIR="/mmfs1/gscratch/cmt/hxd/RealSpace_ExactDiagonalization"
JULIA_BIN="/mmfs1/gscratch/cmt/hxd/opt/julia-1.12.6/bin/julia"
JULIA_DEPOT="/mmfs1/gscratch/cmt/hxd/julia_depot"
ACCOUNT="cmt"
PARTITION="ckpt-g2"
WALLTIME="04:00:00"
MAIL_USER="hxd.phys@outlook.com"

# One generated Slurm script is one independently resumable data point. The
# generated submission helper calls `sbatch` in a loop; `sbatch` returns
# immediately, so all data jobs are queued asynchronously.

# These are numerator values x. The Julia outputs and all plot x-axes use the
# physical hopping t'' = x/(2+2sqrt(2)).
SWEEP_NUMERATORS=(
  -3.0 -2.9 -2.8 -2.7 -2.6 -2.5 -2.4 -2.3 -2.2 -2.1
  -2.0 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1
  -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1
   0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9
   1.0  1.1  1.2  1.3  1.4  1.5
)
PHASES=(AHC FCI CDW)
SWEEP_GEOMETRIES=(3x5 3x6 3x7 4x6)
GAP_GEOMETRIES=(3x4 3x5 3x6 3x7 4x6)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENERATED_DIR="${SCRIPT_DIR}/generated"
HYAK_LOG_DIR="${REPO_DIR}/phase_exploration/hpc/logs"
mkdir -p "${GENERATED_DIR}"
rm -f "${GENERATED_DIR}"/*.sbatch \
      "${GENERATED_DIR}/submit_all.sh" \
      "${GENERATED_DIR}/data_jobs.txt"
DATA_MANIFEST="${GENERATED_DIR}/data_jobs.txt"
: > "${DATA_MANIFEST}"

value_tag() {
  local tag="$1"
  tag="${tag//-/m}"
  tag="${tag//./p}"
  printf '%s' "${tag}"
}

resource_for() {
  local task="$1"
  local geometry="$2"
  MODE="matrix"
  CPUS=1
  MEM_GB=128
  WALLTIME="04:00:00"

  case "${geometry}" in
    3x4) MEM_GB=12 ;;
    3x5) MEM_GB=36 ;;
    3x6) MEM_GB=72 ;;
    3x7) MEM_GB=96 ;;
    4x6) MEM_GB=256 ;;
    *) echo "Unknown geometry ${geometry}" >&2; exit 2 ;;
  esac

  if [[ "${geometry}" == "3x4" ]]; then
    MODE="matrixfree"
    CPUS=8
    WALLTIME="01:00:00"
  elif [[ "${geometry}" == "3x5" ]]; then
    MODE="matrixfree"
    CPUS=12
    WALLTIME="02:00:00"
  elif [[ "${geometry}" == "3x6" ]]; then
    MODE="matrixfree"
    CPUS=72
    WALLTIME="04:00:00"
  elif [[ "${geometry}" == "3x7" ]]; then
    MODE="matrixfree"
    CPUS=144
    WALLTIME="04:00:00"
  elif [[ "${geometry}" == "4x6" ]]; then
    MODE="matrixfree"
    CPUS=216
    WALLTIME="08:00:00"
  fi
}

write_header() {
  local file="$1"
  local job_name="$2"
  cat > "${file}" <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --export=ALL
#SBATCH --job-name=${job_name}
#SBATCH --time=${WALLTIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM_GB}G
#SBATCH --output=${HYAK_LOG_DIR}/${job_name}_%j.out
#SBATCH --error=${HYAK_LOG_DIR}/${job_name}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${MAIL_USER}

set -euo pipefail
export JULIA_PROJECT="${REPO_DIR}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT}"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

WORKERS=\$((SLURM_CPUS_PER_TASK - 1))
JULIA_ARGS=(--project="${REPO_DIR}" --startup-file=no)
if (( WORKERS > 0 )); then
  JULIA_ARGS+=(-p "\${WORKERS}")
fi
EOF
}

for geometry in "${SWEEP_GEOMETRIES[@]}"; do
  resource_for sweep "${geometry}"
  for x in "${SWEEP_NUMERATORS[@]}"; do
    xtag="$(value_tag "${x}")"
    job="${GENERATED_DIR}/sweep_${geometry}_x_${xtag}.sbatch"
    write_header "${job}" "tpp_sw_${geometry}_${xtag}"
    cat >> "${job}" <<EOF
srun "${JULIA_BIN}" "\${JULIA_ARGS[@]}" \
  "${REPO_DIR}/phase_exploration/bin/run_sweep_point.jl" \
  --geometry "${geometry}" --x "${x}" --task all --mode "${MODE}" \
  --nev 10 --dense-resolution 101
EOF
    printf '%s\n' "$(basename "${job}")" >> "${DATA_MANIFEST}"
  done
done

for geometry in "${SWEEP_GEOMETRIES[@]}"; do
  resource_for diagnostics "${geometry}"
  for phase in "${PHASES[@]}"; do
    phase_lower="${phase,,}"
    job="${GENERATED_DIR}/diagnostics_${geometry}_${phase_lower}.sbatch"
    write_header "${job}" "tpp_dx_${geometry}_${phase_lower}"
    cat >> "${job}" <<EOF
srun "${JULIA_BIN}" "\${JULIA_ARGS[@]}" \
  "${REPO_DIR}/phase_exploration/bin/run_diagnostic_point.jl" \
  --phase "${phase}" --geometry "${geometry}" --mode "${MODE}" \
  --observables flow,pump,spatial_es,pes --zero-nev 10 --flow-nev 3 \
  --flow-steps 25 --pump-steps 17 --pes-na 2
EOF
    printf '%s\n' "$(basename "${job}")" >> "${DATA_MANIFEST}"
  done
done

for geometry in "${GAP_GEOMETRIES[@]}"; do
  resource_for charge_gap "${geometry}"
  for phase in "${PHASES[@]}"; do
    phase_lower="${phase,,}"
    job="${GENERATED_DIR}/charge_gap_${geometry}_${phase_lower}.sbatch"
    write_header "${job}" "tpp_cg_${geometry}_${phase_lower}"
    cat >> "${job}" <<EOF
srun "${JULIA_BIN}" "\${JULIA_ARGS[@]}" \
  "${REPO_DIR}/phase_exploration/bin/run_charge_gap_point.jl" \
  --phase "${phase}" --geometry "${geometry}" --mode "${MODE}" --nev 2
EOF
    printf '%s\n' "$(basename "${job}")" >> "${DATA_MANIFEST}"
  done
done

cat > "${GENERATED_DIR}/plot_results.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --export=ALL
#SBATCH --job-name=tpp_plot
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=${HYAK_LOG_DIR}/tpp_plot_%j.out
#SBATCH --error=${HYAK_LOG_DIR}/tpp_plot_%j.err
set -euo pipefail
export JULIA_PROJECT="${REPO_DIR}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT}"
srun "${JULIA_BIN}" --project="${REPO_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/plot_results.jl" --kind all
EOF

cat > "${GENERATED_DIR}/submit_all.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$(dirname "${SCRIPT_DIR}")/logs"
MANIFEST="${SCRIPT_DIR}/data_jobs.txt"
[[ -s "${MANIFEST}" ]] || { echo "Missing or empty ${MANIFEST}" >&2; exit 2; }

mapfile -t job_names < "${MANIFEST}"
job_ids=()
for job_name in "${job_names[@]}"; do
  [[ -n "${job_name}" ]] || continue
  job_path="${SCRIPT_DIR}/${job_name}"
  [[ -f "${job_path}" ]] || { echo "Missing generated job ${job_path}" >&2; exit 2; }
  submission="$(sbatch --parsable "${job_path}")"
  job_id="${submission%%;*}"
  job_ids+=("${job_id}")
  printf 'submitted %-55s -> %s\n' "${job_name}" "${job_id}"
done

if ((${#job_ids[@]} == 0)); then
  echo "No data jobs were submitted." >&2
  exit 2
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
submission_log="${SCRIPT_DIR}/submission_${timestamp}.csv"
{
  echo "job_file,job_id"
  paste -d, "${MANIFEST}" <(printf '%s\n' "${job_ids[@]}")
} > "${submission_log}"

dependency="$(IFS=:; echo "${job_ids[*]}")"
plot_submission="$(sbatch --parsable --dependency="afterok:${dependency}" "${SCRIPT_DIR}/plot_results.sbatch")"
plot_id="${plot_submission%%;*}"
printf 'submitted %-55s -> %s (afterok all data jobs)\n' "plot_results.sbatch" "${plot_id}"
printf 'Queued %d independent data jobs asynchronously.\n' "${#job_ids[@]}"
printf 'Submission record: %s\n' "${submission_log}"
EOF

chmod +x "${GENERATED_DIR}/submit_all.sh"
data_job_count="$(wc -l < "${DATA_MANIFEST}")"
echo "Generated ${data_job_count} independent data jobs plus one dependent plot job in ${GENERATED_DIR}"
echo "Review the collection, then submit ALL jobs asynchronously with:"
echo "  ${GENERATED_DIR}/submit_all.sh"
