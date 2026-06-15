function [static_sinkage_mm, body_weight_N, ground_pressure_Pa] = calc_static_sinkage(...
    mass_kg,...
    gravity_mps2,...
    contact_width_mm,...
    contact_length_mm,...
    kc,...
    kphi,...
    n)
%CALC_STATIC_SINKAGE Bekker式から機体自重による静的沈下量を求める。
%
% 入力
%   mass_kg           : 機体質量 [kg]
%   gravity_mps2      : 重力加速度 [m/s^2]
%   contact_width_mm  : 接地幅 [mm]
%   contact_length_mm : 接地長 [mm]
%   kc, kphi, n       : Bekker圧力沈下式のパラメータ
%
% 出力（返り値）
%   static_sinkage_mm : 静的沈下量 [mm]
%   body_weight_N     : 機体重量 [N]
%   ground_pressure_Pa: 接地圧 [Pa]
%
% 注意
%   この関数は静的平衡のみを扱う。
%   動的な沈下や滑りは考慮しない。

% --------------------------------------------------
% 1. 単位変換
% --------------------------------------------------
% 入力検証
if contact_width_mm <= 0 || contact_length_mm <= 0
    error('contact_width_mm と contact_length_mm は正の値で指定してください（mm）。');
end

contact_width_m = contact_width_mm / 1000;
contact_length_m = contact_length_mm / 1000;
contact_area_m2 = contact_width_m * contact_length_m;

% --------------------------------------------------
% 2. 機体重量と接地圧を計算
% --------------------------------------------------
if mass_kg <= 0 || gravity_mps2 <= 0
    error('mass_kg と gravity_mps2 は正の値を指定してください。');
end

body_weight_N = mass_kg * gravity_mps2;
if contact_area_m2 <= eps
    error('接地面積が 0 または極小です。contact_width_mm と contact_length_mm を確認してください。');
end

ground_pressure_Pa = body_weight_N / contact_area_m2;

% --------------------------------------------------
% 3. Bekker式を逆算して静的沈下量を求める
% --------------------------------------------------
% Bekker式:
%   p = (kc / b + kphi) * z^n
% ここでは b に接地幅を代表値として使う。
% z は m 単位で求め、最後に mm に戻す。
% Bekker 式の分母が 0 にならないように保護
representative_width_m = contact_width_m;
denom = kc / representative_width_m + kphi;
if denom <= eps
    warning('Bekker 式の係数が不適切です。分母を eps にクランプします。');
    denom = max(denom, eps);
end

if n <= 0
    warning('沈下指数 n は正である必要があります。n=1 にフォールバックします。');
    n = 1.0;
end

static_sinkage_m = (ground_pressure_Pa / denom)^(1 / n);
static_sinkage_mm = static_sinkage_m * 1000;

end