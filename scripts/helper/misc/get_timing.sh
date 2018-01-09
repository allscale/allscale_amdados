#!/bin/bash
echo "C++ simulation"

s1=$(date -d"2017-12-25 18:21:09" +%s) # field_Nx121_Ny121_Nt2160_time00000.txt
s2=$(date -d"2017-12-25 18:31:33" +%s) # field_Nx121_Ny121_Nt2160_time02159.txt
echo $((s2-s1))

s1=$(date -d"2018-01-02 19:46:49" +%s) # field_Nx209_Ny187_Nt3546_time00000.txt
s2=$(date -d"2018-01-02 20:21:10" +%s) # field_Nx209_Ny187_Nt3546_time03545.txt
echo $((s2-s1))

s1=$(date -d"2017-12-25 19:26:19" +%s) # field_Nx253_Ny275_Nt4734_time00000.txt
s2=$(date -d"2017-12-25 20:46:41" +%s) # field_Nx253_Ny275_Nt4734_time04733.txt
echo $((s2-s1))

s1=$(date -d"2018-01-02 23:28:01" +%s) # field_Nx407_Ny341_Nt6714_time00000.txt
s2=$(date -d"2018-01-03 03:06:25" +%s) # field_Nx407_Ny341_Nt6714_time06713.txt
echo $((s2-s1))

s1=$(date -d"2017-12-26 04:21:47" +%s) # field_Nx473_Ny451_Nt8298_time00000.txt
s2=$(date -d"2017-12-26 11:07:19" +%s) # field_Nx473_Ny451_Nt8298_time08297.txt
echo $((s2-s1))

s1=$(date -d"2017-12-31 02:16:06" +%s) # field_Nx913_Ny979_Nt17217_time00000.txt
s2=$(date -d"2018-01-02 12:33:38" +%s) # field_Nx913_Ny979_Nt17217_time17216.txt
echo $((s2-s1))

echo
echo
echo "Python simulation"

s1=$(date -d"2017-12-25 18:17:58" +%s) # true-field_Nx121_Ny121_Nt2160_time00000.txt
s2=$(date -d"2017-12-25 18:21:08" +%s) # true-field_Nx121_Ny121_Nt2160_time02159.txt
echo $((s2-s1))

s1=$(date -d"2018-01-02 19:28:04" +%s) # true-field_Nx209_Ny187_Nt3546_time00000.txt
s2=$(date -d"2018-01-02 19:46:44" +%s) # true-field_Nx209_Ny187_Nt3546_time03545.txt
echo $((s2-s1))

s1=$(date -d"2017-12-25 18:31:36" +%s) # true-field_Nx253_Ny275_Nt4734_time00000.txt
s2=$(date -d"2017-12-25 19:26:09" +%s) # true-field_Nx253_Ny275_Nt4734_time04733.txt
echo $((s2-s1))

s1=$(date -d"2018-01-02 20:21:18" +%s) # true-field_Nx407_Ny341_Nt6714_time00000.txt
s2=$(date -d"2018-01-02 23:27:40" +%s) # true-field_Nx407_Ny341_Nt6714_time06713.txt
echo $((s2-s1))

s1=$(date -d"2017-12-25 20:46:55" +%s) # true-field_Nx473_Ny451_Nt8298_time00000.txt
s2=$(date -d"2017-12-26 04:21:08" +%s) # true-field_Nx473_Ny451_Nt8298_time08297.txt
echo $((s2-s1))

s1=$(date -d"2017-12-26 11:07:59" +%s) # true-field_Nx913_Ny979_Nt17217_time00000.txt
s2=$(date -d"2017-12-31 02:11:44" +%s) # true-field_Nx913_Ny979_Nt17217_time17216.txt
echo $((s2-s1))

