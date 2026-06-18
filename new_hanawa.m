%% ==========================
%  初期設定
%% ==========================

% 関連する定数は構造体でまとめて管理する。
limits.Th1_max = deg2rad(225);
limits.Th1_min = deg2rad(-45);
limits.Th2_max = deg2rad(135);
limits.Th2_min = deg2rad(0);

trajectory.y_start = 125;      % 開始位置
trajectory.y_end = 175;        % 終了位置
trajectory.x_fixed = 0;        % x座標固定

% --------------------------------------------------
% Bekker圧力沈下式の条件
% --------------------------------------------------
% ここで求める静的沈下量は、1回だけ計算して全軌道点で共通利用する。
% 接地寸法は mm で定義し、関数内部で m に変換する。
soil.body_mass = 10;         % 機体質量 [kg]
soil.g = 9.81;               % 重力加速度 [m/s^2]
soil.contact_width = 100;    % 接地幅 [mm]
soil.contact_length = 100;   % 接地長 [mm]
soil.kc = 0;                 % Bekker式の凝集項
soil.kphi = 1500e3;          % Bekker式の摩擦項
soil.n = 1.0;                % 沈下指数

[static_sinkage, body_weight, ~] =  calc_static_sinkage( ...
    soil.body_mass, ...
    soil.g, ...
    soil.contact_width, ...
    soil.contact_length, ...
    soil.kc, ...
    soil.kphi, ...
    soil.n);

fprintf('Static sinkage : %.3f mm\n', static_sinkage);
fprintf('Body weight    : %.3f N\n', body_weight);

geometry.L_sum = 250;

%% ==========================
%  配列の準備
%% ==========================

candidate.L1 = zeros(geometry.L_sum-1,1);
candidate.L2 = zeros(geometry.L_sum-1,1);

evaluation.R_all = zeros(geometry.L_sum-1,1);

trajectory.point_num = abs(trajectory.y_start-trajectory.y_end)+1;
trajectory.x = zeros(trajectory.point_num,1);
trajectory.y = zeros(trajectory.point_num,1);
trajectory.sinkage_step = zeros(trajectory.point_num,1);

position.x = zeros(trajectory.point_num,1);
position.y = zeros(trajectory.point_num,1);
position.dx = zeros(trajectory.point_num,1);
position.dy = zeros(trajectory.point_num,1);
position.beta = zeros(trajectory.point_num,1);
position.gamma = zeros(trajectory.point_num,1);

force.alpha_y = zeros(trajectory.point_num,1);
force.alpha_x = zeros(trajectory.point_num,1);
force.Fy = zeros(trajectory.point_num,1);
force.Fx = zeros(trajectory.point_num,1);
force.force_angle = zeros(trajectory.point_num,1);

ellipse.a = zeros(trajectory.point_num,1);
ellipse.b = zeros(trajectory.point_num,1);
ellipse.u1 = zeros(2,trajectory.point_num);
ellipse.u2 = zeros(2,trajectory.point_num);
ellipse.angle = zeros(trajectory.point_num,1);
ellipse.Th = zeros(trajectory.point_num,1);
ellipse.r = zeros(trajectory.point_num,1);

evaluation = struct();

joint = struct();
joint.current.Th1 = 0;
joint.current.Th2 = 0;
joint.next.Th1 = 0;
joint.next.Th2 = 0;

center = struct();
center.current.L2_center = 0;
center.next.L2_center = 0;

sinkage = struct();
sinkage.current = 0;
sinkage.next = 0;

%% ==========================
%  軌道生成
%% ==========================

for i = 1:trajectory.point_num
    
    trajectory.x(i) = trajectory.x_fixed;
    trajectory.y(i) = trajectory.y_start + i - 1;
    trajectory.sinkage_step(i) = i - 1;

end

%% ==========================
%  全リンク長を探索（片方が0はない）
%% ==========================

for c = 1:geometry.L_sum-1
    
    candidate.L1(c) = c;
    candidate.L2(c) = geometry.L_sum - c;
    evaluation.flag = 1;

    %% 到達判定　　...は改行という意味

    for g = 1:trajectory.point_num
        evaluation.flag = judge(...
            trajectory.x(g),...
            trajectory.y(g),...
            candidate.L1(c),...
            candidate.L2(c));

        if evaluation.flag ~= 1
            break
        end

    end

% 到達判定を通過した候補だけを評価する。
    %% 到達可能なら評価値計算
    evaluation.flag=1;
    if evaluation.flag == 1

        evaluation.R = 0;

        for z = 1:trajectory.point_num

            %--------------------------
            % 逆運動学
            %--------------------------
            [joint.current.Th1,joint.current.Th2] = inverse_kinematics(...
                trajectory.x(z),...
                trajectory.y(z),...
                candidate.L1(c),...
                candidate.L2(c));

            if any(~isfinite([joint.current.Th1, joint.current.Th2]))
                evaluation.flag = -3;
                break
            end

            % 関節角制限違反なら、このリンク長候補を失格にする。
            if angle_judge(joint.current.Th1, limits.Th1_min, limits.Th1_max) ~= 1 || ...
                    angle_judge(joint.current.Th2, limits.Th2_min, limits.Th2_max) ~= 1
                evaluation.flag = -2;
                break
            end

            %--------------------------
            % 地盤への貫入量
            %--------------------------

            % Bekker式から求めた静的沈下量と、軌道に沿った追加沈下量を足す。
            % sinkage.current は、目標軌道が開始位置より何 mm 深いかを表す。
            sinkage.current = static_sinkage + (trajectory.y(z) - trajectory.y_start);

            % 地中部分の中心位置

            % Th1, Th2 はラジアンなので、cos() には pi を使う。
            center.current.L2_center = candidate.L2(c) ...
                - sinkage.current * cos(pi - joint.current.Th1 - joint.current.Th2) / 2;

            %--------------------------
            % 順運動学
            %--------------------------

            [position.x(z),position.y(z)] = ...
                forward_kinematics(...
                joint.current.Th1,...
                joint.current.Th2,...
                candidate.L1(c),...
                center.current.L2_center);

            %--------------------------
            % 移動方向
            %--------------------------

            if z == trajectory.point_num

                position.dx(z) = 0;
                position.dy(z) = 0;

            else

                % 次の点を計算

                

            %--------------------------
            % RFT（砂地盤では、足部の姿勢や移動方向によって地盤反力が変化するため、RFT（Resistive Force Theory）を用いて砂から受ける反力を推定する。
            %--------------------------
                [joint.next.Th1,joint.next.Th2] = ...
                    inverse_kinematics(...
                    trajectory.x(z+1),...
                    trajectory.y(z+1),...
                    candidate.L1(c),...
                    candidate.L2(c));

                if any(~isfinite([joint.next.Th1, joint.next.Th2]))
                    evaluation.flag = -3;
                    break
                end

                sinkage.next = ...
                    static_sinkage + (trajectory.y(z+1) - trajectory.y_start);

                center.next.L2_center = ...
                    candidate.L2(c) - sinkage.next ...
                    * cos(pi - joint.next.Th1 - joint.next.Th2) / 2;

                [position.next.x,position.next.y] = ...
                    forward_kinematics(...
                    joint.next.Th1,...
                    joint.next.Th2,...
                    candidate.L1(c),...
                    center.next.L2_center);

                position.dx(z) = position.next.x-position.x(z);
                position.dy(z) = position.next.y-position.y(z);

            end  %ここまででx、yの位置、x、yの移動量

            %--------------------------
            % β・γ(アームの角度と進行方向の向きの決定をしている)
            %--------------------------

            position.beta(z) = joint.current.Th1 + joint.current.Th2;
            position.gamma(z) = atan2(...
                position.dy(z),...
                position.dx(z));

            force.alpha_y(z) = 0.055*sin(-2*position.beta(z)+position.gamma(z))+0.206+0.358*sin(position.gamma(z))+0.169*cos(2*position.beta(z))+0.212*sin(2*position.beta(z)+position.gamma(z));

            force.alpha_x(z) = -0.124*cos(2*position.beta(z)+position.gamma(z))+0.253*cos(position.gamma(z))+0.007*cos(-2*position.beta(z)+position.gamma(z))+0.088*sin(2*position.beta(z));

            force.Fy(z) = 0.191 * force.alpha_y(z) * trajectory.sinkage_step(z);
            force.Fx(z) = 0.191 * force.alpha_x(z) * trajectory.sinkage_step(z);

            force.force_angle(z) = atan2(force.Fy(z),force.Fx(z));  %反力の向き

            %--------------------------
            % 浮き上がり判定
            %--------------------------
            % 鉛直方向のみを判定する。
            % Fy が body_weight を超えた時点で、底面反力は 0 未満になれないため
            % 機体が浮き始めると判断する。
            if force.Fy(z) > body_weight
                % このリンク長候補は失格。残りの軌道計算は打ち切る。
                evaluation.flag = -5;
                break
            end

            %--------------------------
            % 操作力楕円体
            %--------------------------

            [ellipse.a(z),ellipse.b(z),ellipse.u1(:,z),ellipse.u2(:,z)] = ...
                ellipse_parameters(...
                candidate.L1(c),...
                center.current.L2_center,...
                joint.current.Th1,...
                joint.current.Th2);

            ellipse.angle(z) = ...
                atan2(ellipse.u2(2,z),ellipse.u2(1,z));

            ellipse.Th(z) = force.force_angle(z) ...
                + ellipse.angle(z);  %力の作用線と操作力楕円体の主軸との角度

            ellipse.r(z) = ellipse_r(ellipse.a(z),ellipse.b(z),ellipse.Th(z));

            % 評価値更新
            
            evaluation.R = evaluation.R + ellipse.r(z);
            
        end

        % 1つでも条件違反があれば、この c 候補の残り処理を飛ばして次へ進む。
        if evaluation.flag ~= 1
            evaluation.R = -Inf;
            evaluation.R_all(c) = evaluation.R;
            continue
        end

    end

    evaluation.R_all(c) = evaluation.R;

end

%% ==========================
% 最適リンク長
%% ==========================

[M,I] = max(evaluation.R_all);

if isinf(M) && M < 0
    fprintf("L1,L2の候補がありません\n");
else
    L1 = candidate.L1(I);
    L2 = candidate.L2(I);
    fprintf("L1 = %d\n", L1);
    fprintf("L2 = %d\n", L2);
end
%% ==========================
% function
%% ==========================

function [Th1,Th2]=inverse_kinematics(x,y,L1,L2)
    % inverse_kinematics 2リンク平面ロボットの逆運動学（肘下解）
    % 入力:
    %   x, y  - エンドエフェクタ目標位置（L1,L2 と同じ単位）
    %   L1,L2 - リンク長
    % 出力:
    %   Th1, Th2 - 関節角（ラジアン）
    % 注意:
    %   acos の引数は [-1,1] に収める必要があります（到達不能点で複素数が発生します）。
    reach = (x^2+y^2+L1^2-L2^2)/(2*L1*sqrt(x^2+y^2));
    if ~isreal(reach) || abs(reach) > 1
        Th1 = NaN;
        Th2 = NaN;
        return
    end

    Th1=atan2(y,x)-acos(reach);
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
