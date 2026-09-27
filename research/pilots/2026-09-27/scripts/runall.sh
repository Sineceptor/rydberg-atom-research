run(){ python scan.py "$@" >> results.jsonl 2>>err.log; }
( run 0.6 11.39 2.8 4.5; run 0.3 11.39 0.5 3.0 ) &
( run 0.4 11.39 0.8 3.5; run 0.7 40 0.5 6.0 ) &
( run 0.3 40 0.3 4.0; run 0.5 150 0.5 8.0 ) &
( run 0.5 40 0.5 5.0; run 0.4 100 0.5 6.0 ) &
wait
