for NAME in $(screen -ls | grep [0-9][0-9][0-9][0-9])
do
    echo $NAME
    screen -S $NAME -X quit
done

killall MATLAB