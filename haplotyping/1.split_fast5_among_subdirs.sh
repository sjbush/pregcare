n=0; for f in *.fast5; do d="subdir$((n++ / 3))"; mkdir -p "$d"; mv -- "$f" "$d/$f"; done

