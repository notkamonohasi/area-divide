target="main.cpp"

echo
g++ -O3 -o temp_main ${target}
echo "compile done"
echo

if test -e temp_main; then
    start_time=$(gdate +'%H%M%S.%3N')
    ./temp_main <in.txt >out.txt
    end_time=$(gdate +'%H%M%S.%3N')
    calcTime=$(echo "scale=2; ($end_time - $start_time)" | bc | xargs printf "%.2f\n")
    echo
    echo calculation Time : ${calcTime} s
    echo
    rm temp_main
fi
