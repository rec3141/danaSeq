#!/usr/bin/env bash
# Record the promoted digest of an image in deploy/<image>.json on the branch,
# which omc-pickup reads to pin each cluster to an exact image.
#
#   record-deploy.sh <image> <digest> [base_digest] [built]
#
# Several builds can finish close together (a base build chains every leaf,
# and pushes to one leaf can overlap), so each attempt starts from the current
# branch head and rewrites the file rather than rebasing onto another build's
# version of it. A build never replaces the record of a newer commit.
set -euo pipefail
image=$1 digest=$2 base=${3:-} built=${4:-}
sha=${GITHUB_SHA:?} ref=${GITHUB_REF_NAME:?}
file="deploy/${image}.json"

git config user.name  "github-actions[bot]"
git config user.email "github-actions[bot]@users.noreply.github.com"

# The ancestry test below needs history; CI checkouts are shallow.
if [ "$(git rev-parse --is-shallow-repository)" = "true" ]; then
    git fetch -q --unshallow origin "$ref" || true
fi

for i in 1 2 3 4 5; do
    git fetch -q origin "$ref"
    git checkout -q -B "$ref" FETCH_HEAD
    recorded=$(jq -r '.commit // empty' "$file" 2>/dev/null || true)
    if [ -n "$recorded" ] && [ "$recorded" != "$sha" ]; then
        git fetch -q origin "$recorded" 2>/dev/null || true
        if git merge-base --is-ancestor "$sha" "$recorded" 2>/dev/null; then
            echo "$file already records $recorded, which is newer than $sha; leaving it"
            exit 0
        fi
    fi
    mkdir -p deploy
    jq -n --arg image "$image" --arg digest "$digest" --arg commit "$sha" \
          --arg base "$base" --arg built "$built" \
          '{image: $image, digest: $digest, commit: $commit}
           + (if $base  != "" then {base_digest: $base} else {} end)
           + (if $built != "" then {built: $built} else {} end)' > "$file"
    git add "$file"
    if git diff --cached --quiet; then echo "digest unchanged"; exit 0; fi
    git commit -q -m "deploy: ${image} ${sha} [skip ci]"
    git push -q origin "HEAD:${ref}" && { echo "recorded ${image} ${digest}"; exit 0; }
    sleep $((RANDOM % 5 + 2))
done
echo "could not push the deploy record after 5 attempts" >&2
exit 1
