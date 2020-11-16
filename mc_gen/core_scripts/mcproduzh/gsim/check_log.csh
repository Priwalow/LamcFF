#!/bin/csh -f

# This is a script for checking the contents in log files.
#
# Usage:
# 1. Put the log files named "*.log" at the directory named "log".
# 2. Exec "./check_log.csh".
# 3. Check the characters shown.
#    Different numbers for $1 at "grep" enable you to check several values.
#    (If you are satisfied here, the task of this script is finished.)
# 4. Modify "???" in the file "check_log.csh"
#    referring the result at the step 3.
# 5. Exec "./check_log.csh".
# 6. Check the log files listed in the output "list_failed.txt".
#    If there are no problems, the output "list_failed.txt" has no characters.
#
# Contact: Y. Horii (yhorii@epx.phys.tohoku.ac.jp)
# Date: 24.3.2010

set list_all = "list_all.txt"
set list_ok = "list_ok.txt"
set list_failed = "list_failed.txt"
rm $list_all
rm $list_ok
rm $list_failed

ls log/*log >>! $list_all

foreach f(log/*.log)
 set grp = `grep "\#output events" $f | awk '{print $1}'` # $1 can be changed.
                                                          # Try below.
# set grp = `grep -a "\#output events" $f | awk '{print $NF}'` # $NF can be changed.
# set grp = `grep -a "\#input events" $f | awk '{print $NF}'` # $NF can be changed.
 echo $grp
 if ($grp == "???") then    # ??? should be modified at the step 4.
  echo $f >>! $list_ok
 endif
end
echo "" >>! $list_ok

comm -23 $list_all $list_ok >>! $list_failed

rm $list_all
rm $list_ok
