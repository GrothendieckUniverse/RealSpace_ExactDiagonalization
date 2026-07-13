#!/usr/bin/env bash

set -euo pipefail

# ==========================================================================
# User-editable Hyak / Klone paths and allocation
# ==========================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# By default, infer the repository root from this generator's location. This
# keeps generated jobs valid when the checkout lives under, for example,
# /gscratch/cmt/hxd/repo instead of a hard-coded /mmfs1 path. Export REPO_DIR
# before running this script only if an explicit override is needed.
REPO_DIR="${REPO_DIR:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
JULIA_BIN="${JULIA_BIN:-/mmfs1/gscratch/cmt/hxd/opt/julia-1.12.6/bin/julia}"
JULIA_DEPOT="${JULIA_DEPOT:-/mmfs1/gscratch/cmt/hxd/julia_depot}"
JULIA_PROJECT_DIR="${JULIA_PROJECT_DIR:-${JULIA_DEPOT}/environments/v1.12}"
ACCOUNT="cmt"
PARTITION="ckpt-g2"
WALLTIME="04:00:00"
MAIL_USER="hxd.phys@outlook.com"
REPO_REVISION="$(git -C "${REPO_DIR}" rev-parse --short=8 HEAD 2>/dev/null || printf 'working')"
SETUP_JOB_NAME="tpp_env_${REPO_REVISION}"

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
SWEEP_GEOMETRIES=(3x5 3x6 3x7)
GAP_GEOMETRIES=(3x4 3x5 3x6 3x7)

GENERATED_DIR="${SCRIPT_DIR}/generated"
HYAK_LOG_DIR="${REPO_DIR}/phase_exploration/hpc/logs"
mkdir -p "${GENERATED_DIR}" "${HYAK_LOG_DIR}"
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

result_value_tag() {
  local formatted
  printf -v formatted '%.4f' "$1"
  value_tag "${formatted}"
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
    3x6) MEM_GB=48 ;;
    3x7) MEM_GB=72 ;;
    4x6) MEM_GB=240 ;;
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
    CPUS=36
    WALLTIME="04:00:00"
  elif [[ "${geometry}" == "3x7" ]]; then
    MODE="matrixfree"
    CPUS=72
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
#SBATCH --ntasks-per-node=${CPUS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_GB}G
#SBATCH --chdir=${REPO_DIR}
#SBATCH --output=${HYAK_LOG_DIR}/${job_name}_%j.out
#SBATCH --error=${HYAK_LOG_DIR}/${job_name}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${MAIL_USER}

set -Eeuo pipefail
trap 'status=\$?; printf "ERROR: exit=%d line=%d command=%s\\n" "\${status}" "\${LINENO}" "\${BASH_COMMAND}" >&2; exit "\${status}"' ERR
export JULIA_PROJECT="${JULIA_PROJECT_DIR}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT}"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export SLURM_EXPORT_ENV=ALL
export SRUN_EXPORT_ENV=ALL
export PHASE_EXPLORATION_REPO="${REPO_DIR}"
DONE_DIR="${REPO_DIR}/phase_exploration/hpc/completed"
DONE_FILE="\${DONE_DIR}/${job_name}.done"

echo "Starting Slurm job \${SLURM_JOB_ID:-unknown} on \$(hostname) at \$(date --iso-8601=seconds)"
echo "Repository: ${REPO_DIR}"
echo "Julia project: ${JULIA_PROJECT_DIR}"
echo "Julia: ${JULIA_BIN}"
echo "Resources: tasks=\${SLURM_NTASKS:-unknown} memory=${MEM_GB}G"
[[ -d "${REPO_DIR}" ]] || { echo "Missing repository: ${REPO_DIR}" >&2; exit 2; }
[[ -f "${REPO_DIR}/Project.toml" ]] || { echo "Missing Project.toml under ${REPO_DIR}" >&2; exit 2; }
[[ -f "${JULIA_PROJECT_DIR}/Project.toml" ]] || { echo "Missing shared Julia environment: ${JULIA_PROJECT_DIR}" >&2; exit 2; }
[[ -x "${JULIA_BIN}" ]] || { echo "Julia is not executable: ${JULIA_BIN}" >&2; exit 2; }
"${JULIA_BIN}" --version

mark_complete() {
  mkdir -p "\${DONE_DIR}"
  printf 'job_id=%s\\ncompleted_at=%s\\n' \
    "\${SLURM_JOB_ID:-unknown}" "\$(date --iso-8601=seconds)" > "\${DONE_FILE}.tmp"
  mv "\${DONE_FILE}.tmp" "\${DONE_FILE}"
  echo "Completion marker: \${DONE_FILE}"
}
EOF
}

cat > "${GENERATED_DIR}/setup_environment.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --export=ALL
#SBATCH --job-name=${SETUP_JOB_NAME}
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --chdir=${REPO_DIR}
#SBATCH --output=${HYAK_LOG_DIR}/${SETUP_JOB_NAME}_%j.out
#SBATCH --error=${HYAK_LOG_DIR}/${SETUP_JOB_NAME}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${MAIL_USER}

set -Eeuo pipefail
trap 'status=\$?; printf "ERROR: exit=%d line=%d command=%s\n" "\${status}" "\${LINENO}" "\${BASH_COMMAND}" >&2; exit "\${status}"' ERR
export JULIA_PROJECT="${JULIA_PROJECT_DIR}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT}"
export PHASE_EXPLORATION_REPO="${REPO_DIR}"
DONE_DIR="${REPO_DIR}/phase_exploration/hpc/completed"
DONE_FILE="\${DONE_DIR}/${SETUP_JOB_NAME}.done"

echo "Preparing ${JULIA_PROJECT_DIR} for repository ${REPO_DIR}"
[[ -x "${JULIA_BIN}" ]] || { echo "Julia is not executable: ${JULIA_BIN}" >&2; exit 2; }
"${JULIA_BIN}" --project="${JULIA_PROJECT_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/prepare_hpc_environment.jl"
mkdir -p "\${DONE_DIR}"
printf 'job_id=%s\ncompleted_at=%s\nrevision=%s\n' \
  "\${SLURM_JOB_ID:-unknown}" "\$(date --iso-8601=seconds)" "${REPO_REVISION}" > "\${DONE_FILE}"
echo "Environment marker: \${DONE_FILE}"
EOF

for geometry in "${SWEEP_GEOMETRIES[@]}"; do
  resource_for sweep "${geometry}"
  for x in "${SWEEP_NUMERATORS[@]}"; do
    xtag="$(value_tag "${x}")"
    result_xtag="$(result_value_tag "${x}")"
    job="${GENERATED_DIR}/sweep_${geometry}_x_${xtag}.sbatch"
    write_header "${job}" "tpp_sw_${geometry}_${xtag}"
    cat >> "${job}" <<EOF
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/sweep/${geometry}/x_${result_xtag}/spectrum.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/sweep/${geometry}/x_${result_xtag}/structure_allowed.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/sweep/${geometry}/x_${result_xtag}/structure_dense.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/sweep/${geometry}/x_${result_xtag}/structure_metrics.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/sweep/${geometry}/x_${result_xtag}/run_metadata.csv
"${JULIA_BIN}" --project="${JULIA_PROJECT_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/run_slurm_job.jl" \
  "${REPO_DIR}/phase_exploration/bin/run_sweep_point.jl" \
  --geometry "${geometry}" --x "${x}" --task all --mode "${MODE}" \
  --nev 10 --dense-resolution 101
mark_complete
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
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/zero_flux_spectrum.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/spectrum_flow.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/charge_pump.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/spatial_entanglement_spectrum.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/particle_entanglement_spectrum.csv
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/diagnostics/${phase}/${geometry}/summary.csv
"${JULIA_BIN}" --project="${JULIA_PROJECT_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/run_slurm_job.jl" \
  "${REPO_DIR}/phase_exploration/bin/run_diagnostic_point.jl" \
  --phase "${phase}" --geometry "${geometry}" --mode "${MODE}" \
  --observables flow,pump,spatial_es,pes --zero-nev 10 --flow-nev 3 \
  --flow-steps 25 --pump-steps 17 --pes-na 2
mark_complete
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
# PHASE_STUDY_REQUIRED_OUTPUT=${REPO_DIR}/phase_exploration/results/charge_gap/${phase}/${geometry}/charge_gap.csv
"${JULIA_BIN}" --project="${JULIA_PROJECT_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/run_slurm_job.jl" \
  "${REPO_DIR}/phase_exploration/bin/run_charge_gap_point.jl" \
  --phase "${phase}" --geometry "${geometry}" --mode "${MODE}" --nev 2
mark_complete
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
#SBATCH --chdir=${REPO_DIR}
#SBATCH --output=${HYAK_LOG_DIR}/tpp_plot_%j.out
#SBATCH --error=${HYAK_LOG_DIR}/tpp_plot_%j.err
set -Eeuo pipefail
trap 'status=\$?; printf "ERROR: exit=%d line=%d command=%s\\n" "\${status}" "\${LINENO}" "\${BASH_COMMAND}" >&2; exit "\${status}"' ERR
export JULIA_PROJECT="${JULIA_PROJECT_DIR}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT}"
export PHASE_EXPLORATION_REPO="${REPO_DIR}"
DONE_DIR="${REPO_DIR}/phase_exploration/hpc/completed"
DONE_FILE="\${DONE_DIR}/tpp_plot.done"
echo "Starting plot job \${SLURM_JOB_ID:-unknown} on \$(hostname) at \$(date --iso-8601=seconds)"
[[ -d "${REPO_DIR}" ]] || { echo "Missing repository: ${REPO_DIR}" >&2; exit 2; }
[[ -x "${JULIA_BIN}" ]] || { echo "Julia is not executable: ${JULIA_BIN}" >&2; exit 2; }
"${JULIA_BIN}" --version
"${JULIA_BIN}" --project="${REPO_DIR}" --startup-file=no \
  "${REPO_DIR}/phase_exploration/bin/plot_results.jl" --kind all
mkdir -p "\${DONE_DIR}"
printf 'job_id=%s\\ncompleted_at=%s\\n' \
  "\${SLURM_JOB_ID:-unknown}" "\$(date --iso-8601=seconds)" > "\${DONE_FILE}"
echo "Completion marker: \${DONE_FILE}"
EOF

cat > "${GENERATED_DIR}/submit_all.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HPC_DIR="$(dirname "${SCRIPT_DIR}")"
DONE_DIR="${HPC_DIR}/completed"
mkdir -p "${HPC_DIR}/logs" "${DONE_DIR}"
MANIFEST="${SCRIPT_DIR}/data_jobs.txt"
[[ -s "${MANIFEST}" ]] || { echo "Missing or empty ${MANIFEST}" >&2; exit 2; }

declare -A active_ids_by_name=()
while IFS='|' read -r active_name active_id; do
  [[ -n "${active_name}" && -n "${active_id}" ]] || continue
  if [[ -n "${active_ids_by_name[${active_name}]:-}" ]]; then
    active_ids_by_name["${active_name}"]+=":${active_id}"
  else
    active_ids_by_name["${active_name}"]="${active_id}"
  fi
done < <(squeue --noheader --user "${USER}" --format='%j|%A' 2>/dev/null || true)

job_slurm_name() {
  local job_path="$1"
  sed -n 's/^#SBATCH --job-name=//p' "${job_path}"
}

job_outputs_complete() {
  local job_path="$1"
  local required_output
  local found=0
  while IFS= read -r required_output; do
    [[ -n "${required_output}" ]] || continue
    found=1
    [[ -s "${required_output}" ]] || return 1
  done < <(sed -n 's/^# PHASE_STUDY_REQUIRED_OUTPUT=//p' "${job_path}")
  (( found == 1 ))
}

declare -A dependency_seen=()
dependency_ids=()
add_dependency_ids() {
  local id
  local ids="$1"
  local split_ids=()
  IFS=':' read -r -a split_ids <<< "${ids}"
  for id in "${split_ids[@]}"; do
    [[ -n "${id}" ]] || continue
    if [[ -z "${dependency_seen[${id}]:-}" ]]; then
      dependency_seen["${id}"]=1
      dependency_ids+=("${id}")
    fi
  done
}

mapfile -t job_names < "${MANIFEST}"
job_ids=()
submission_records=()
skipped_active=0
skipped_complete=0

setup_job_path="${SCRIPT_DIR}/setup_environment.sbatch"
[[ -f "${setup_job_path}" ]] || { echo "Missing ${setup_job_path}" >&2; exit 2; }
setup_slurm_name="$(job_slurm_name "${setup_job_path}")"
[[ -n "${setup_slurm_name}" ]] || { echo "Missing Slurm job name in ${setup_job_path}" >&2; exit 2; }
setup_marker="${DONE_DIR}/${setup_slurm_name}.done"
setup_dependency=""

if [[ -n "${active_ids_by_name[${setup_slurm_name}]:-}" ]]; then
  setup_dependency="${active_ids_by_name[${setup_slurm_name}]}"
  add_dependency_ids "${setup_dependency}"
  submission_records+=("setup_environment.sbatch,skipped_active,${setup_dependency//:/|}")
  printf 'skipped   %-55s -> environment setup already active as %s\n' \
    "setup_environment.sbatch" "${setup_dependency}"
elif [[ -f "${setup_marker}" ]]; then
  submission_records+=("setup_environment.sbatch,skipped_complete,environment_marker")
  printf 'skipped   %-55s -> environment marker already exists\n' \
    "setup_environment.sbatch"
else
  setup_submission="$(sbatch --parsable "${setup_job_path}")"
  setup_dependency="${setup_submission%%;*}"
  active_ids_by_name["${setup_slurm_name}"]="${setup_dependency}"
  add_dependency_ids "${setup_dependency}"
  submission_records+=("setup_environment.sbatch,submitted,${setup_dependency}")
  printf 'submitted %-55s -> %s\n' "setup_environment.sbatch" "${setup_dependency}"
fi

for job_name in "${job_names[@]}"; do
  [[ -n "${job_name}" ]] || continue
  job_path="${SCRIPT_DIR}/${job_name}"
  [[ -f "${job_path}" ]] || { echo "Missing generated job ${job_path}" >&2; exit 2; }
  slurm_name="$(job_slurm_name "${job_path}")"
  [[ -n "${slurm_name}" ]] || { echo "Missing Slurm job name in ${job_path}" >&2; exit 2; }

  if [[ -n "${active_ids_by_name[${slurm_name}]:-}" ]]; then
    active_ids="${active_ids_by_name[${slurm_name}]}"
    add_dependency_ids "${active_ids}"
    submission_records+=("${job_name},skipped_active,${active_ids//:/|}")
    printf 'skipped   %-55s -> already queued/running as %s\n' "${job_name}" "${active_ids}"
    ((skipped_active += 1))
    continue
  fi

  if job_outputs_complete "${job_path}"; then
    printf 'detected_from=result_files\ncompleted_at=%s\n' \
      "$(date --iso-8601=seconds)" > "${DONE_DIR}/${slurm_name}.done"
    submission_records+=("${job_name},skipped_complete,result_files")
    printf 'skipped   %-55s -> required result files already exist\n' "${job_name}"
    ((skipped_complete += 1))
    continue
  fi

  sbatch_args=(--parsable)
  if [[ -n "${setup_dependency}" ]]; then
    sbatch_args+=(--dependency="afterok:${setup_dependency}")
  fi
  submission="$(sbatch "${sbatch_args[@]}" "${job_path}")"
  job_id="${submission%%;*}"
  job_ids+=("${job_id}")
  add_dependency_ids "${job_id}"
  active_ids_by_name["${slurm_name}"]="${job_id}"
  submission_records+=("${job_name},submitted,${job_id}")
  printf 'submitted %-55s -> %s\n' "${job_name}" "${job_id}"
done

timestamp="$(date +%Y%m%d_%H%M%S)"
submission_log="${SCRIPT_DIR}/submission_${timestamp}.csv"
{
  echo "job_file,status,job_id_or_reason"
  printf '%s\n' "${submission_records[@]}"
} > "${submission_log}"

if [[ -n "${active_ids_by_name[tpp_plot]:-}" ]]; then
  printf 'skipped   %-55s -> plot job already queued/running as %s\n' \
    "plot_results.sbatch" "${active_ids_by_name[tpp_plot]}"
elif ((${#dependency_ids[@]} > 0)); then
  dependency="$(IFS=:; echo "${dependency_ids[*]}")"
  plot_submission="$(sbatch --parsable --dependency="afterok:${dependency}" "${SCRIPT_DIR}/plot_results.sbatch")"
  plot_id="${plot_submission%%;*}"
  printf 'submitted %-55s -> %s (afterok active and newly submitted data jobs)\n' \
    "plot_results.sbatch" "${plot_id}"
elif [[ -f "${DONE_DIR}/tpp_plot.done" ]]; then
  printf 'skipped   %-55s -> completion marker already exists\n' "plot_results.sbatch"
else
  plot_submission="$(sbatch --parsable "${SCRIPT_DIR}/plot_results.sbatch")"
  plot_id="${plot_submission%%;*}"
  printf 'submitted %-55s -> %s (all data results already complete)\n' \
    "plot_results.sbatch" "${plot_id}"
fi

printf 'Submitted %d new data jobs; skipped %d active and %d complete jobs.\n' \
  "${#job_ids[@]}" "${skipped_active}" "${skipped_complete}"
printf 'Submission record: %s\n' "${submission_log}"
EOF

chmod +x "${GENERATED_DIR}/submit_all.sh"
data_job_count="$(wc -l < "${DATA_MANIFEST}")"
echo "Generated one environment job, ${data_job_count} independent data jobs, and one dependent plot job in ${GENERATED_DIR}"
echo "Review the collection, then submit ALL jobs asynchronously with:"
echo "  ${GENERATED_DIR}/submit_all.sh"
