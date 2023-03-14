for N in 20 35; do
    mkdir Nc$N
    cp Nc50/input.jl Nc$N/
    sed -i "s/Nc = 50/Nc = $N/g" Nc$N/input.jl
done
grep 'Nc = ' */input.jl
echo '------------------'
grep 'σM =' */input.jl

