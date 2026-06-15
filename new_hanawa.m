%% ==========================
%  初期設定
%% ==========================

% 関節角度制限
global Th1_max Th1_min Th2_max Th2_min y1 y2 x1

Th1_max = 225;
Th1_min = -45;
Th2_max = 135;
Th2_min = 0;

% 軌道条件
y1 = 125;      % 開始位置
y2 = 225;      % 終了位置
x1 = 0;        % x座標固定

% --------------------------------------------------
% Bekker圧力沈下式の条件
% --------------------------------------------------
% ここで求める静的沈下量は、1回だけ計算して全軌道点で共通利用する。
% 接地寸法は mm で定義し、関数内部で m に変換する。
body_mass = 0.01;         % 機体質量 [kg]
g = 9.81;               % 重力加速度 [m/s^2]
contact_width = 100;    % 接地幅 [mm]
contact_length = 100;   % 接地長 [mm]
kc = 0;                 % Bekker式の凝集項
kphi = 1500e3;          % Bekker式の摩擦項
n = 1.0;                % 沈下指数

[static_sinkage, body_weight, ~] =  calc_static_sinkage( ...
    body_mass, ...
    g, ...
    contact_width, ...
    contact_length, ...
    kc, ...
    kphi, ...
    n);

fprintf('Static sinkage : %.3f mm\n', static_sinkage);
fprintf('Body weight    : %.3f N\n', body_weight);

% リンク長の総和
L_sum = 250;

%% ==========================
%  配列の準備
%% ==========================

con = zeros(L_sum-1,2);     % リンク長候補
R_all = zeros(L_sum-1,1);   % 評価値

point_num = abs(y1-y2)+1;

tra = zeros(point_num,3);   %1mmずつの沈下の回数
pos = zeros(point_num,6);

%% ==========================
%  軌道生成
%% ==========================

for i = 1:point_num

    tra(i,1) = x1;
    tra(i,2) = y1 + i - 1;
    tra(i,3) = i - 1;

end

%% ==========================
%  全リンク長を探索（片方が0はない）
%% ==========================

for c = 1:L_sum-1

    L1 = c;
    L2 = L_sum - c;

    con(c,1) = L1;   %con リンク長候補の配列
    con(c,2) = L2;

    %% 到達判定　　...は改行という意味

    for g = 1:point_num

        flag = judge(...
            tra(g,1),...
            tra(g,2),...
            L1,...
            L2);

        if flag ~= 1
            break
        end

    end
%50回の判定処理は終了している。一度でも判定が失敗になったらそのアームの組み合わせ除外外される
    %% 到達可能なら評価値計算

    if flag == 1

        R = 0;

        for z = 1:point_num

            %--------------------------
            % 逆運動学
            %--------------------------

            [Th1,Th2] = inverse_kinematics(...
                tra(z,1),...
                tra(z,2),...
                L1,...
                L2);

            %--------------------------
            % 地盤への貫入量
            %--------------------------

            % Bekker式から求めた静的沈下量と、軌道に沿った追加沈下量を足す。
            % motion_sinkage は、目標軌道が開始位置より何 mm 深いかを表す。
            motion_sinkage = tra(z,2) - y1;
            sinkage = static_sinkage + motion_sinkage;

            % 地中部分の中心位置

            % Th1, Th2 はラジアンなので、cos() には pi を使う。
            L2_center = L2 ...
                - sinkage * cos(pi - Th1 - Th2) / 2;

            %--------------------------
            % 順運動学
            %--------------------------

            [pos(z,1),pos(z,2)] = ...
                forward_kinematics(...
                Th1,...
                Th2,...
                L1,...
                L2_center);

            %--------------------------
            % 移動方向
            %--------------------------

            if z == point_num

                pos(z,3) = 0;
                pos(z,4) = 0;

            else

                % 次の点を計算

                [Th1_next,Th2_next] = ...
                    inverse_kinematics(...
                    tra(z+1,1),...
                    tra(z+1,2),...
                    L1,...
                    L2);

                motion_sinkage_next = ...
                    tra(z+1,2) - y1;

                sinkage_next = ...
                    static_sinkage + motion_sinkage_next;

                L2_center_next = ...
                    L2 - sinkage_next ...
                    * cos(pi - Th1_next - Th2_next) / 2;

                [x_next,y_next] = ...
                    forward_kinematics(...
                    Th1_next,...
                    Th2_next,...
                    L1,...
                    L2_center_next);

                pos(z,3) = x_next-pos(z,1);
                pos(z,4) = y_next-pos(z,2);

            end  %ここまででx、yの位置、x、yの移動量アームと体の向きの決定をしている

            %--------------------------
            % β・γ
            %--------------------------

            beta = Th1 + Th2;
            gamma = atan2(...
                pos(z,4),...
                pos(z,3));

            pos(z,5) = beta;
            pos(z,6) = gamma;

            %--------------------------
            % RFT（砂地盤では、足部の姿勢や移動方向によって地盤反力が変化するため、RFT（Resistive Force Theory）を用いて砂から受ける反力を推定する。
            %--------------------------

            alpha_y = ...
                0.055*sin(-2*beta+gamma) ...
                +0.206 ...
                +0.358*sin(gamma) ...
                +0.169*cos(2*beta) ...
                +0.212*sin(2*beta+gamma);

            alpha_x = ...
                -0.124*cos(2*beta+gamma) ...
                +0.253*cos(gamma) ...
                +0.007*cos(-2*beta+gamma) ...
                +0.088*sin(2*beta);

            Fy = 0.191 * alpha_y * tra(z,3);
            Fx = 0.191 * alpha_x * tra(z,3);

            force_angle = atan2(Fy,Fx);

            %--------------------------
            % 浮き上がり判定
            %--------------------------
            % 鉛直方向のみを判定する。
            % Fy が body_weight を超えた時点で、底面反力は 0 未満になれないため
            % 機体が浮き始めると判断する。
            if Fy > body_weight
                flag = -5;
                break
            end

            %--------------------------
            % 操作力楕円体
            %--------------------------

            [a,b,u1,u2] = ...
                ellipse_parameters(...
                L1,...
                L2_center,...
                Th1,...
                Th2);

            ellipse_angle = ...
                atan2(u2(2),u2(1));

            Th = force_angle ...
                + ellipse_angle;

            r = ellipse_r(a,b,Th);

            % 評価値更新

            R = R + r;

        end

    else

        R = flag;

    end

    R_all(c) = R;

end

%% ==========================
% 最適リンク長
%% ==========================

[M,I] = max(R_all);

L1 = con(I,1);
L2 = con(I,2);

disp(['L1 : ',num2str(L1)])
disp(['L2 : ',num2str(L2)])


function[Th1,Th2]=inverse_kinematics(x,y,L1,L2)
    % inverse_kinematics 2リンク平面ロボットの逆運動学（肘下解）
    % 入力:
    %   x, y  - エンドエフェクタ目標位置（L1,L2 と同じ単位）
    %   L1,L2 - リンク長
    % 出力:
    %   Th1, Th2 - 関節角（ラジアン）
    % 注意:
    %   acos の引数は [-1,1] に収める必要があります（到達不能点で複素数が発生します）。
    Th1=atan2(y,x)-acos((x^2+y^2+L1^2-L2^2)/(2*L1*sqrt(x^2+y^2)));
    Th2=-Th1+atan2((y-L1*sin(Th1)),(x-L1*cos(Th1)));
    Th1=normalize_angle(Th1);
    Th2=normalize_angle(Th2);
end

function [x, y] = forward_kinematics(Th1, Th2, L1, L2)
    % forward_kinematics 2リンク平面ロボットの順運動学
    % 入力:
    %   Th1, Th2 - 関節角（ラジアン）
    %   L1, L2   - リンク長
    % 出力:
    %   x, y - エンドエフェクタ位置
    x = L1*cos(Th1) + L2*cos(Th1 + Th2);
    y = L1*sin(Th1) + L2*sin(Th1 + Th2);
end

function[a,b,u1,u2]=ellipse_parameters(L1,L2,Th1,Th2)
    % ellipse_parameters 操作力（または速度）楕円体のパラメータを計算
    % 入力:
    %   L1,L2,Th1,Th2 - リンク長と関節角
    % 出力:
    %   a,b   - 楕円の長短半径（SVD の特異値の逆数）
    %   u1,u2 - 楕円の主軸方向（列ベクトル）
    J = [
        -L1*sin(Th1)-L2*sin(Th1+Th2),-L2*sin(Th1+Th2);
         L1*cos(Th1)+L2*cos(Th1+Th2), L2*cos(Th1+Th2);
    ];
    [U,S,~]=svd(J);
    a=1/S(1,1);
    b=1/S(2,2);
    u1=U(:,1);
    u2=U(:,2);
end

function [r]=ellipse_r(a,b,Th)
    % ellipse_r 楕円上の方向角 Th における半径を返す
    % 入力:
    %   a,b - 楕円の半径（長軸 a, 短軸 b）
    %   Th  - 角度（ラジアン）
    % 出力:
    %   r - その方向の距離
    x=(a*b)/sqrt(b^2+(a*tan(Th))^2);
    y=(a*b*tan(Th))/sqrt(b^2+(a*tan(Th))^2);
    r=sqrt(x^2+y^2);
end

function normalized_angle=normalize_angle(angle)
    % normalize_angle 角度を範囲 [-pi/2, 3*pi/2) に正規化する
    % 入力:
    %   angle - 任意の角度（ラジアン）
    % 出力:
    %   normalized_angle - 正規化後の角度
    if angle>3*pi/2||angle<-pi/2
        normalized_angle=mod(angle+pi/2,2*pi)-pi/2;
    else
        normalized_angle=angle;
    end
end
