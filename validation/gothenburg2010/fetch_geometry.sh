#!/usr/bin/env sh
set -eu

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
data_directory="${script_directory}/data"
archive_directory="${data_directory}/archives"

mkdir -p "${archive_directory}"
mkdir -p "${data_directory}/KVLCC2"
mkdir -p "${data_directory}/KCS"

download_and_verify() {
    url=$1
    archive=$2
    checksum=$3

    checksum_line="${checksum}  ${archive}"
    if [ -f "${archive}" ] && \
        printf '%s\n' "${checksum_line}" | sha256sum --check --status
    then
        return
    fi

    partial_archive="${archive}.part"
    curl --fail --location --retry 3 --output "${partial_archive}" "${url}"
    partial_checksum_line="${checksum}  ${partial_archive}"
    printf '%s\n' "${partial_checksum_line}" | sha256sum --check --status
    mv "${partial_archive}" "${archive}"
}

base_url="https://www.nmri.go.jp/archives/institutes/fluid_performance_evaluation/cfd_rd/cfdws05/gothenburg2000/data"
kvlcc2_archive="${archive_directory}/kvlcc2_surfacemesh.zip"
kcs_archive="${archive_directory}/kcs_surfacemesh.zip"

download_and_verify \
    "${base_url}/KVLCC2_data/kvlcc2_surfacemesh.ZIP" \
    "${kvlcc2_archive}" \
    "5d2b8023ec0efdb0e8c3c979a0e0d68ead56cd8c77e31a4da3500a06800b2d2a"
download_and_verify \
    "${base_url}/KCS_data/kcs_surfacemesh.ZIP" \
    "${kcs_archive}" \
    "cccefd3aca810b5ec6691f494fe254b1c5ef9e9329758900d4d64683206ed71a"

unzip -oq "${kvlcc2_archive}" -d "${data_directory}/KVLCC2"
unzip -oq "${kcs_archive}" -d "${data_directory}/KCS"

printf '%s\n' "Workshop hull surfaces are available in ${data_directory}."
