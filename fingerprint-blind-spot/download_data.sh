#!/usr/bin/env bash
# Fetch the MoleculeNet CSVs used in this analysis and check SHA-256.
# Source: DeepChem public S3 bucket. These files are not stored in the repo.
set -euo pipefail
cd "$(dirname "$0")"
BASE=https://deepchemdata.s3-us-west-1.amazonaws.com/datasets

fetch_and_check() {
  local file="$1" expected="$2"
  echo "downloading $file"
  curl -fL --retry 3 -o "$file" "$BASE/$file"
  local got
  got="$(python3 -c "import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" "$file")"
  if [ "$got" != "$expected" ]; then
    echo "SHA-256 mismatch for $file" >&2
    echo "  expected $expected" >&2
    echo "  got      $got" >&2
    exit 1
  fi
  echo "ok $file"
}

fetch_and_check delaney-processed.csv 8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189
fetch_and_check SAMPL.csv             ab5895d914ee87cb563bd7b9611e869527bba45bec6b014d34dc495a0f9dcb72
fetch_and_check Lipophilicity.csv     aed41590cb30609d51d8e08ad3ff06495a76e80e211358801f596b10da69bacd
fetch_and_check bace.csv              f3fb9ce90bada3e2bd6148b0df13f8f8145a357bf87df0dd5b391ede974fc737
fetch_and_check BBBP.csv              d07a38487aeac5cee5508413e468043ef3097451d2a112701c2d60be9ec6b662
echo "all five files match the hashes in the README"
