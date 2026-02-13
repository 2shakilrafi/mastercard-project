
#!/usr/bin/env bash

shopt -s nullglob

for f in *; do
  new="$(echo "$f" | tr '[:upper:]' '[:lower:]' | sed 's/ /-/g')"

  # Skip if name is unchanged
  [[ "$f" == "$new" ]] && continue

  # Skip if target already exists
  if [[ -e "$new" ]]; then
    echo "Skipping '$f' → '$new' (target exists)"
    continue
  fi

  mv -- "$f" "$new"
done
