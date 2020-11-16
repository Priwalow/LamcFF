#!/usr/bin/expect -f
set timeout -1
spawn ./bbsview
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "L\r"
expect "*\r"
send -- "Q\r"
expect eof 
