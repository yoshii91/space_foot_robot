function flag = judge(x, y, L1, L2)
%JUDGE 2リンク平面アームが目標点に到達可能かを判定する。
%
% 入力
%   x, y : 目標点座標 [mm]
%   L1   : 第1リンク長 [mm]
%   L2   : 第2リンク長 [mm]
%
% 出力
%   flag : 1 なら到達可能、-1 なら到達不能

distance = hypot(x, y);     %sqrt(x^2 + y^2)
reach_max = L1 + L2;
reach_min = abs(L1 - L2);
tolerance = 1e-9 * max([reach_max, 1]);

if distance <= reach_max + tolerance && distance >= reach_min - tolerance
    flag = 1;
else
    flag = -1;
end

end