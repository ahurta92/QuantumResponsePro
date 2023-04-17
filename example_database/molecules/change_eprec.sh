#/bin/bash
#!/usr/bin/env bash

for i in *.mol; do
    [ -f "$i" ] || break
    sed -i 's/eprec 1e-7/eprec 1e-6/' $i
done
