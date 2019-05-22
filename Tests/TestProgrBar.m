% T0024: Test case for Progress bar
% Johan Gustafsson, 2019-05-21
% See VerificationMatrix.csv for how test cases map to code

%The expected output is a progress bar called "B1 - sub", that progresses
%at different speeds. It should progress to 60%, then halt for a while,
%then go faster to 75%, then go slower to 90%, then halt for a while, then
%progress to 100. No other bars should be shown.

b1 = ProgrBar('B1 - sub');
b2 = ProgrBar('B2', b1.GetSubContext(0.6));
b3 = ProgrBar('B3', b1.GetSubContext(0.3));
b4 = ProgrBar('B4', b3.GetSubContext(0.5));
b5 = ProgrBar('B5', b3.GetSubContext(0.5));
%b5 = ProgrBar('B5', b3, 0.6);%replace the above line with this one to test that the sum of fractions should not go above 1

for i = 1:100
    pause(0.05)
    b2.Progress(i/100);
end
b2.Done();

pause(1)

for i = 1:100
    pause(0.01)
    b4.Progress(i/100);
end
b4.Done();

for i = 1:100
    pause(0.05)
    b5.Progress(i/100);
end
b5.Done();

b3.Done();

pause(1)

for i = 1:100
    pause(0.02)
    b1.Progress(i/100);
end


b1.Done();
b1.Done(); %check that one extra Done doesn't mess things up

%test silence
bs = ProgrBar('BS - This should not be shown', ProgrBar.GetSilentContext());
for i = 1:100
    pause(0.02)
    bs.Progress(i/100);
end

bs.Done();


