mkdir Nc75
mkdir Nc100
cp Nc50/input.jl Nc75/
cp Nc50/input.jl Nc100/
sed -i 's/Nc = 50/Nc = 75/g' Nc75/input.jl
sed -i 's/Nc = 50/Nc = 100/g' Nc100/input.jl
grep 'Nc = ' */input.jl
echo '------------------'
grep 'σM =' */input.jl

