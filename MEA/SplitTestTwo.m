function [a,b,c] = SplitTestTwo(filename1,filename2)

TempOne = importtest(filename1);
TempTwo = importtest(filename2);
%%
%WT
%A1 A2 B1 B2 C1 C2 D1 D2
%00 01 06 07 12 13 18 19
%01 02 07 08 13 14 19 20

a.filename = TempOne.filename;

a.Well(1:2) = TempOne.Well(1:2);
a.Well(3:4) = TempOne.Well(7:8);
a.Well(5:6) = TempOne.Well(13:14);
a.Well(7:8) = TempOne.Well(19:20);
a.Well(9:10) = TempTwo.Well(1:2);
a.Well(11:12) = TempTwo.Well(7:8);
a.Well(13:14) = TempTwo.Well(13:14);
a.Well(15:16) = TempTwo.Well(19:20);

a.syncIndex(1:2) = TempOne.syncIndex(1:2);
a.syncIndex(3:4) = TempOne.syncIndex(7:8);
a.syncIndex(5:6) = TempOne.syncIndex(13:14);
a.syncIndex(7:8) = TempOne.syncIndex(19:20);
a.syncIndex(9:10) = TempTwo.syncIndex(1:2);
a.syncIndex(11:12) = TempTwo.syncIndex(7:8);
a.syncIndex(13:14) = TempTwo.syncIndex(13:14);
a.syncIndex(15:16) = TempTwo.syncIndex(19:20);

%%
%KO
%A3 A4 B3 B4 C3 C4 D3 D4
%02 03 08 09 14 15 20 21
%03 04 09 10 15 16 21 22

b.filename = TempOne.filename;

b.Well(1:2) = TempOne.Well(3:4);
b.Well(3:4) = TempOne.Well(9:10);
b.Well(5:6) = TempOne.Well(15:16);
b.Well(7:8) = TempOne.Well(21:22);
b.Well(9:10) = TempTwo.Well(3:4);
b.Well(11:12) = TempTwo.Well(9:10);
b.Well(13:14) = TempTwo.Well(15:16);
b.Well(15:16) = TempTwo.Well(21:22);

b.syncIndex(1:2) = TempOne.syncIndex(3:4);
b.syncIndex(3:4) = TempOne.syncIndex(9:10);
b.syncIndex(5:6) = TempOne.syncIndex(15:16);
b.syncIndex(7:8) = TempOne.syncIndex(21:22);
b.syncIndex(9:10) = TempTwo.syncIndex(3:4);
b.syncIndex(11:12) = TempTwo.syncIndex(9:10);
b.syncIndex(13:14) = TempTwo.syncIndex(15:16);
b.syncIndex(15:16) = TempTwo.syncIndex(21:22);


%%
%MIX
%A5 A6 B5 B6 C5 C6 D5 D6
%04 05 10 11 16 17 22 23
%05 06 11 12 17 18 23 24

c.filename = TempOne.filename;

c.Well(1:2) = TempOne.Well(5:6);
c.Well(3:4) = TempOne.Well(11:12);
c.Well(5:6) = TempOne.Well(17:18);
c.Well(7:8) = TempOne.Well(23:24);
c.Well(9:10) = TempTwo.Well(5:6);
c.Well(11:12) = TempTwo.Well(11:12);
c.Well(13:14) = TempTwo.Well(17:18);
c.Well(15:16) = TempTwo.Well(23:24);

c.syncIndex(1:2) = TempOne.syncIndex(5:6);
c.syncIndex(3:4) = TempOne.syncIndex(11:12);
c.syncIndex(5:6) = TempOne.syncIndex(17:18);
c.syncIndex(7:8) = TempOne.syncIndex(23:24);
c.syncIndex(9:10) = TempTwo.syncIndex(5:6);
c.syncIndex(11:12) = TempTwo.syncIndex(11:12);
c.syncIndex(13:14) = TempTwo.syncIndex(17:18);
c.syncIndex(15:16) = TempTwo.syncIndex(23:24);