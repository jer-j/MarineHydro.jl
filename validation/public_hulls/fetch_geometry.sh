#!/usr/bin/env sh
set -eu

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
data_directory="${script_directory}/data"

mkdir -p "${data_directory}/DTMB5415"
mkdir -p "${data_directory}/DTC"

download_and_verify() {
    url=$1
    destination=$2
    checksum=$3

    checksum_line="${checksum}  ${destination}"
    if [ -f "${destination}" ] && \
        printf '%s\n' "${checksum_line}" | sha256sum --check --status
    then
        return
    fi

    partial_destination="${destination}.part"
    curl --fail --location --retry 3 --output "${partial_destination}" "${url}"
    partial_checksum_line="${checksum}  ${partial_destination}"
    printf '%s\n' "${partial_checksum_line}" | sha256sum --check --status
    mv "${partial_destination}" "${destination}"
}

nmri_base="https://www.nmri.go.jp/archives/institutes/fluid_performance_evaluation/cfd_rd/cfdws05/gothenburg2000/data/5415"
dtmb_surface="${data_directory}/DTMB5415/5415_static.net"
download_and_verify \
    "${nmri_base}/5415_static.net" \
    "${dtmb_surface}" \
    "6cda8ea29d4efecab814737d32d487916ed7ba5cbbc85242e77a44ba36612a34"

dtc_archive="${data_directory}/DTC/DTC-scaled.stl.gz"
download_and_verify \
    "https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-dev/e11dbc66d29b34031e9cc6335ac20ae2d1487af2/tutorials/resources/geometry/DTC-scaled.stl.gz" \
    "${dtc_archive}" \
    "61a60ffca0aae2ffb7703cdaf1bab4de7dfeb719ccc3ea8ae9eadd1932c097a0"
gzip -dkf "${dtc_archive}"

"${script_directory}/../gothenburg2010/fetch_geometry.sh"
printf '%s\n' "Public hull surfaces are available in ${data_directory}."
