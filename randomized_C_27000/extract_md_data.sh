#|/usr/bin/usr
#set -xe

get_files(){
    n=1
    filename=$1
    ext="${filename##*.}"
    name="${filename%.${ext}}"

    str=""
    while [ -f "${name}_${n}.${ext}" ]; do
	str="${str} ${name}_${n}.${ext}"
        n=$((n+1))
    done
    str="${str} $filename"

    echo $str
}

get_data(){
    strings=( time energy energy_exp temperature )

    infile=$1

    prefix=$2

    first_line="#"
    all_pastes=""

    for i in ${!strings[@]}; do
	grep -oP "${strings[$i]}=-?\d+\.\d+" $infile | cut -d '=' -f 2 > temp_$i
	all_pastes="${all_pastes} temp_$i"
	first_line="${first_line} ${strings[$i]}"
    done

    f="${first_line:2}"
    f="${f// /_}"

    file="${prefix}${f}_${infile}.dat"
    #    echo "${first_line}" > $file
    paste $all_pastes > $file

    echo $file
}

cleanup(){
    rm -f temp_*
#    t=$1
#    rm -f ${t}*
}

file="trajectory_out.xyz"
prefix="md_"

files=$(get_files $file)

out_files=""
n=0
offset=""
out=""
nf=${#files}
for f in $files; do

    out=$(get_data $f $prefix )
    out_files="${out_files} ${out}"
    echo "> Output of $f in $out"

    if [ $n != 0 ]; then
	awk "{ \$1+=${offset}; print }" $out > ${out}_temp
	mv ${out}_temp $out
    fi

    n=$((n+1))

    if [ $n > 0  ]; then
	offset=$(tail -1 $out | cut  -f 1)
    fi
done

echo $out_files
n=0
for f in $out_files; do
    if [ $n == 0 ]; then
	cat $out_files > ${f/.dat/_all.dat}
    fi
    rm $f

    n=$((n+1))
done

cleanup
