#!/usr/bin/expect -f
set timeout -1
spawn bbsview /group/belle/bdata_b/mcprod/dat/e000035/evtgen/mixed/07/all/0127/on_resonance/04/evtgen-mixed-07-all-e000035r000457-b20090127_0910.mdst
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "Q\r"
expect eof 
