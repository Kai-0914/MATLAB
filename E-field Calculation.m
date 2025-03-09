clear all;

%% constant
f1 = 2.4*10^9; 
f2 = 5.8*10^9; 
lightspeed = 3*10^8; 
beta1 = 2*pi*f1/lightspeed; 
beta2 = 2*pi*f2/lightspeed; 

%% coordinate
% general
sourceX = 0; sourceY = 4; sourceZ = 2; % source coordinate
receiverX = 10; receiverY = 3; receiverZ = 2; % p coordinate
PtoQ = receiverX:.01:22; % from p to q

% first order coordinate (top)
I1_topX = 0; I1_topY = 8; I1_topZ = 2; 
reflect1_topY = 6; 
reflect1_topZ = 2;
reflect1_topX = (reflect1_topY - I1_topY) ./ ((receiverY - I1_topY) ./ (receiverX - I1_topX));

% first order coordinate (bottom)
I1_bottomX = 0; I1_bottomY = -4; I1_bottomZ = 2; 
reflect1_bottomY = 0; 
reflect1_bottomZ = 2;
reflect1_bottomX = (reflect1_bottomY - I1_bottomY) ./ ((receiverY - I1_bottomY) ./ (receiverX - I1_bottomX));

% second order coordinate (top)
I2a_topX = 0; I2a_topY = 8; I2a_topZ = 2; % first image
I2b_topX = 0; I2b_topY = -8; I2b_topZ = 2; % second image
reflect2b_topY = 0; 
reflect2b_topZ = 2;
reflect2b_topX = (reflect2b_topY - I2b_topY) ./ ((receiverY - I2b_topY) ./ (receiverX - I2b_topX));
reflect2a_topY = 6; 
reflect2a_topZ = 2;
reflect2a_topX = (reflect2a_topY - I2a_topY) ./ ((I2a_topY - reflect2b_topY) ./ (I2a_topX - reflect2b_topX));
intersect_topX = 4;
intersect_topY = 3.6;
intersect_topZ = 2;

% second order coordinate (bottom)
I2a_bottomX = 0; I2a_bottomY = -4; I2a_bottomZ = 2; % first image
I2b_bottomX = 0; I2b_bottomY = 16; I2b_bottomZ = 2; % second image
reflect2a_bottomY = 0; 
reflect2a_bottomZ = 2;
reflect2a_bottomX = (reflect2a_bottomY - I2b_bottomY) ./ ((receiverY - I2b_bottomY) ./ (receiverX - I2b_bottomX));
reflect2b_bottomY = 6; 
reflect2b_bottomZ = 2;
reflect2b_bottomX = (reflect2b_bottomY - I2a_bottomY) ./ ((I2a_bottomY - reflect2a_bottomY) ./ (I2a_bottomX - reflect2a_bottomX));
intersect_bottomX = 7.7;
intersect_bottomY = 6;
intersect_bottomZ = 2;

% third order coordinate (top)
I3a_topX = 0; I3a_topY = 8; I3a_topZ = 2; % first image
I3b_topX = 0; I3b_topY = -8; I3b_topZ = 2; % second image
I3c_topX = 0; I3c_topY = 14; I3c_topZ = 2; % third image
reflect3c_topY = 6; 
reflect3c_topZ = 2;
reflect3c_topX = (reflect3c_topY - I3c_topY) ./ ((receiverY - I3c_topY) ./ (receiverX - I3c_topX));
reflect3b_topY = 0; 
reflect3b_topZ = 2;
reflect3b_topX = (reflect3b_topY - I3b_topY) ./ ((reflect3c_topY - I3b_topY) ./ (reflect3c_topX - I3b_topX));
reflect3a_topY = 6; 
reflect3a_topZ = 2;
reflect3a_topX = (reflect3a_topY - I3a_topY) ./ ((I3a_topY - reflect3b_topY) ./ (I3a_topX - reflect3b_topX));
intersect3a_topX = 2.19;
intersect3a_topY = 3.78;
intersect3a_topZ = 2;
intersect3b_topX = 5.93;
intersect3b_topY = 3.4;
intersect3b_topZ = 2;

% third order coordinate (bottom)
I3a_bottomX = 0; I3a_bottomY = -4; I3a_bottomZ = 2; % first image
I3b_bottomX = 0; I3b_bottomY = 16; I3b_bottomZ = 2; % second image
I3c_bottomX = 0; I3c_bottomY = -16; I3c_bottomZ = 2; % third image
reflect3c_bottomY = 0; 
reflect3c_bottomZ = 2;
reflect3c_bottomX = (reflect3c_bottomY - I3c_bottomY) ./ ((receiverY - I3c_bottomY) ./ (receiverX - I3c_bottomX));
reflect3b_bottomY = 6; 
reflect3b_bottomZ = 2;
reflect3b_bottomX = (reflect3b_bottomY - I3b_bottomY) ./ ((reflect3c_bottomY - I3b_bottomY) ./ (reflect3c_bottomX - I3b_bottomX));
reflect3a_bottomY = 0; 
reflect3a_bottomZ = 2;
reflect3a_bottomX = (reflect3a_bottomY - I3a_bottomY) ./ ((I3a_bottomY - reflect3b_bottomY) ./ (I3a_bottomX - reflect3b_bottomX));
intersect3a_bottomX = 4;
intersect3a_bottomY = 3.6;
intersect3a_bottomZ = 2;
intersect3b_bottomX = 6.67;
intersect3b_bottomY = 3.33;
intersect3b_bottomZ = 2;

% fourth order coordinate (top)
I4a_topX = 0; I4a_topY = 8; I4a_topZ = 2; % first image
I4b_topX = 0; I4b_topY = -8; I4b_topZ = 2; % second image
I4c_topX = 0; I4c_topY = 14; I4c_topZ = 2; % third image
I4d_topX = 0; I4d_topY = -14; I4d_topZ = 2; % fourth image
reflect4d_topY = 0; 
reflect4d_topZ = 2;
reflect4d_topX = (reflect4d_topY - I4d_topY) ./ ((receiverY - I4d_topY) ./ (receiverX - I4d_topX));
reflect4c_topY = 6; 
reflect4c_topZ = 2;
reflect4c_topX = (reflect4c_topY - I4c_topY) ./ ((I4c_topY - reflect4d_topY) ./ (I4c_topX - reflect4d_topX));
reflect4b_topY = 0; 
reflect4b_topZ = 2;
reflect4b_topX = (reflect4b_topY - I4b_topY) ./ ((I4b_topY - reflect4c_topY) ./ (I4b_topX - reflect4c_topX));
reflect4a_topY = 6; 
reflect4a_topZ = 2;
reflect4a_topX = (reflect4a_topY - I4a_topY) ./ ((I4a_topY - reflect4b_topY) ./ (I4a_topX - reflect4b_topX));
intersect4a_topX = 1.4;
intersect4a_topY = 3.86;
intersect4a_topZ = 2;
intersect4b_topX = 3.88;
intersect4b_topY = 3.61;
intersect4b_topZ = 2;
intersect4c_topX = 6.25;
intersect4c_topY = 3.38;
intersect4c_topZ = 2;

% fourth order coordinate (bottom)
I4a_bottomX = 0; I4a_bottomY = -4; I4a_bottomZ = 2; % first image
I4b_bottomX = 0; I4b_bottomY = 16; I4b_bottomZ = 2; % second image
I4c_bottomX = 0; I4c_bottomY = -16; I4c_bottomZ = 2; % third image
I4d_bottomX = 0; I4d_bottomY = 22; I4d_bottomZ = 2; % fourth image
reflect4d_bottomY = 6; 
reflect4d_bottomZ = 2;
reflect4d_bottomX = (reflect4d_bottomY - I4d_bottomY) ./ ((receiverY - I4d_bottomY) ./ (receiverX - I4d_bottomX));
reflect4c_bottomY = 0; 
reflect4c_bottomZ = 2;
reflect4c_bottomX = (reflect4c_bottomY - I4c_bottomY) ./ ((I4c_bottomY - reflect4d_bottomY) ./ (I4c_bottomX - reflect4d_bottomX));
reflect4b_bottomY = 6; 
reflect4b_bottomZ = 2;
reflect4b_bottomX = (reflect4b_bottomY - I4b_bottomY) ./ ((I4b_bottomY - reflect4c_bottomY) ./ (I4b_bottomX - reflect4c_bottomX));
reflect4a_bottomY = 0; 
reflect4a_bottomZ = 2;
reflect4a_bottomX = (reflect4a_bottomY - I4a_bottomY) ./ ((I4a_bottomY - reflect4b_bottomY) ./ (I4a_bottomX - reflect4b_bottomX));
intersect4a_bottomX = 2.96;
intersect4a_bottomY = 3.7;
intersect4a_bottomZ = 2;
intersect4b_bottomX = 4.8;
intersect4b_bottomY = 3.52;
intersect4b_bottomZ = 2;
intersect4c_bottomX = 7.4;
intersect4c_bottomY = 3.26;
intersect4c_bottomZ = 2;

%% line of sight
distance_direct = sqrt(((PtoQ - sourceX).^2)+((receiverY - sourceY)^2)+((receiverZ - sourceZ)^2)); % distance along x-axis
Ed_1 = (1./distance_direct).*exp(-j * beta1 * distance_direct); % direct ray (line of sight)
Ed_2 = (1./distance_direct).*exp(-j * beta2 * distance_direct);
EddB_1 = 20*log10(abs(Ed_1)); % direct ray alone
EddB_2 = 20*log10(abs(Ed_2)); % direct ray alone

%% first order
% top
distance1_top = sqrt(((PtoQ - I1_topX).^2)+((receiverY - I1_topY)^2)+((receiverZ - I1_topZ)^2)); % distance along x-axis 
A1_top = sqrt(((sourceX - reflect1_topX).^2)+((sourceY - reflect1_topY)^2)+((sourceZ - reflect1_topZ)^2));
B1_top = sqrt(((reflect1_topX - receiverX).^2)+((reflect1_topY - receiverY)^2)+((reflect1_topZ - receiverZ)^2));
C1_top = sqrt(((sourceX - receiverX).^2)+((sourceY - receiverY)^2)+((sourceZ - receiverZ)^2));
theta1_top = acosd((A1_top^2 + B1_top^2 - C1_top^2) / (2 * A1_top * B1_top));
theta1_top_i = theta1_top ./ 2;
theta1_top_t = asind((1/sqrt(3)) * sind(theta1_top_i));
coefficient1_top = (cosd(theta1_top_i) - (sqrt(3) * cosd(theta1_top_t))) / (cosd(theta1_top_i) + (sqrt(3) * cosd(theta1_top_t)));
E1_top_1 = ((1./distance1_top).*exp(-j * beta1 * distance1_top)) * coefficient1_top;
E1_top_2 = ((1./distance1_top).*exp(-j * beta2 * distance1_top)) * coefficient1_top;

% bottom
distance1_bottom = sqrt(((PtoQ - I1_bottomX).^2)+((receiverY - I1_bottomY)^2)+((receiverZ - I1_bottomZ)^2)); % distance along x-axis 
A1_bottom = sqrt(((sourceX - reflect1_bottomX).^2)+((sourceY - reflect1_bottomY)^2)+((sourceZ - reflect1_bottomZ)^2));
B1_bottom = sqrt(((reflect1_bottomX - receiverX).^2)+((reflect1_bottomY - receiverY)^2)+((reflect1_bottomZ - receiverZ)^2));
C1_bottom = sqrt(((sourceX - receiverX).^2)+((sourceY - receiverY)^2)+((sourceZ - receiverZ)^2));
theta1_bottom = acosd((A1_bottom^2 + B1_bottom^2 - C1_bottom^2) / (2 * A1_bottom * B1_bottom));
theta1_bottom_i = theta1_bottom ./ 2;
theta1_bottom_t = asind((1/sqrt(3)) * sind(theta1_bottom_i));
coefficient1_bottom = (cosd(theta1_bottom_i) - (sqrt(3) * cosd(theta1_bottom_t))) / (cosd(theta1_bottom_i) + (sqrt(3) * cosd(theta1_bottom_t)));
E1_bottom_1 = ((1./distance1_bottom).*exp(-j * beta1 * distance1_bottom)) * coefficient1_bottom;
E1_bottom_2 = ((1./distance1_bottom).*exp(-j * beta2 * distance1_bottom)) * coefficient1_bottom;

% total
E1_1 = Ed_1 + E1_top_1 + E1_bottom_1;
E1_2 = Ed_2 + E1_top_2 + E1_bottom_2;
E1dB_1 = 20*log10(abs(E1_1));
E1dB_2 = 20*log10(abs(E1_2));

%% second order
% top
distance2_top = sqrt(((PtoQ - I2b_topX).^2)+((receiverY - I2b_topY)^2)+((receiverZ - I2b_topZ)^2));
A2a_top = sqrt(((sourceX - reflect2a_topX).^2)+((sourceY - reflect2a_topY)^2)+((sourceZ - reflect2a_topZ)^2));
B2a_top = sqrt(((reflect2a_topX - intersect_topX).^2)+((reflect2a_topY - intersect_topY)^2)+((reflect2a_topZ - intersect_topZ)^2));
C2a_top = sqrt(((sourceX - intersect_topX).^2)+((sourceY - intersect_topY)^2)+((sourceZ - intersect_topZ)^2));
A2b_top = sqrt(((intersect_topX - reflect2b_topX).^2)+((intersect_topY - reflect2b_topY)^2)+((intersect_topZ - reflect2b_topZ)^2));
B2b_top = sqrt(((receiverX - reflect2b_topX).^2)+((receiverY - reflect2b_topY)^2)+((receiverZ - reflect2b_topZ)^2));
C2b_top = sqrt(((intersect_topX - receiverX).^2)+((intersect_topY - receiverY)^2)+((intersect_topZ - receiverZ)^2));
theta2a_top = acosd((A2a_top^2 + B2a_top^2 - C2a_top^2) / (2 * A2a_top * B2a_top));
theta2b_top = acosd((A2b_top^2 + B2b_top^2 - C2b_top^2) / (2 * A2b_top * B2b_top));
theta2a_top_i = theta2a_top ./ 2;
theta2b_top_i = theta2b_top ./ 2;
theta2a_top_t = asind((1/sqrt(3)) * sind(theta2a_top_i));
theta2b_top_t = asind((1/sqrt(3)) * sind(theta2b_top_i));
coefficient2a_top = (cosd(theta2a_top_i) - (sqrt(3) * cosd(theta2a_top_t))) / (cosd(theta2a_top_i) + (sqrt(3) * cosd(theta2a_top_t)));
coefficient2b_top = (cosd(theta2b_top_i) - (sqrt(3) * cosd(theta2b_top_t))) / (cosd(theta2b_top_i) + (sqrt(3) * cosd(theta2b_top_t)));
E2_top_1 = ((1./distance2_top).*exp(-j * beta1 * distance2_top)) * (coefficient2a_top * coefficient2b_top);
E2_top_2 = ((1./distance2_top).*exp(-j * beta2 * distance2_top)) * (coefficient2a_top * coefficient2b_top);

% bottom
distance2_bottom = sqrt(((PtoQ - I2b_bottomX).^2)+((receiverY - I2b_bottomY)^2)+((receiverZ - I2b_bottomZ)^2));
A2a_bottom = sqrt(((sourceX - reflect2a_bottomX).^2)+((sourceY - reflect2a_bottomY)^2)+((sourceZ - reflect2a_bottomZ)^2));
B2a_bottom = sqrt(((reflect2a_bottomX - intersect_bottomX).^2)+((reflect2a_bottomY - intersect_bottomY)^2)+((reflect2a_bottomZ - intersect_bottomZ)^2));
C2a_bottom = sqrt(((sourceX - intersect_bottomX).^2)+((sourceY - intersect_bottomY)^2)+((sourceZ - intersect_bottomZ)^2));
A2b_bottom = sqrt(((intersect_bottomX - reflect2b_bottomX).^2)+((intersect_bottomY - reflect2b_bottomY)^2)+((intersect_bottomZ - reflect2b_bottomZ)^2))
B2b_bottom = sqrt(((receiverX - reflect2b_bottomX).^2)+((receiverY - reflect2b_bottomY)^2)+((receiverZ - reflect2b_bottomZ)^2));
C2b_bottom = sqrt(((intersect_bottomX - receiverX).^2)+((intersect_bottomY - receiverY)^2)+((intersect_bottomZ - receiverZ)^2));
theta2a_bottom = acosd((A2a_bottom^2 + B2a_bottom^2 - C2a_bottom^2) / (2 * A2a_bottom * B2a_bottom));
theta2b_bottom = acosd((A2b_bottom^2 + B2b_bottom^2 - C2b_bottom^2) / (2 * A2b_bottom * B2b_bottom));
theta2a_bottom_i = theta2a_bottom ./ 2;
theta2b_bottom_i = theta2b_bottom ./ 2;
theta2a_bottom_t = asind((1/sqrt(3)) * sind(theta2a_bottom_i));
theta2b_bottom_t = asind((1/sqrt(3)) * sind(theta2b_bottom_i));
coefficient2a_bottom = (cosd(theta2a_bottom_i) - (sqrt(3) * cosd(theta2a_bottom_t))) / (cosd(theta2a_bottom_i) + (sqrt(3) * cosd(theta2a_bottom_t)));
coefficient2b_bottom = (cosd(theta2b_bottom_i) - (sqrt(3) * cosd(theta2b_bottom_t))) / (cosd(theta2b_bottom_i) + (sqrt(3) * cosd(theta2b_bottom_t)));
E2_bottom_1 = ((1./distance2_bottom).*exp(-j * beta1 * distance2_bottom)) * (coefficient2a_bottom * coefficient2b_bottom);
E2_bottom_2 = ((1./distance2_bottom).*exp(-j * beta2 * distance2_bottom)) * (coefficient2a_bottom * coefficient2b_bottom);

% total
E2_1 = Ed_1 + E1_top_1 + E1_bottom_1 + E2_top_1 + E2_bottom_1;
E2_2 = Ed_2 + E1_top_2 + E1_bottom_2 + E2_top_2 + E2_bottom_2;
E2dB_1 = 20*log10(abs(E2_1));
E2dB_2 = 20*log10(abs(E2_2));

%% third order
% top
distance3_top = sqrt(((PtoQ - I3c_topX).^2)+((receiverY - I3c_topY)^2)+((receiverZ - I3c_topZ)^2));
A3a_top = sqrt(((sourceX - reflect3a_topX).^2)+((sourceY - reflect3a_topY)^2)+((sourceZ - reflect3a_topZ)^2));
B3a_top = sqrt(((reflect3a_topX - intersect3a_topX).^2)+((reflect3a_topY - intersect3a_topY)^2)+((reflect3a_topZ - intersect3a_topZ)^2));
C3a_top = sqrt(((sourceX - intersect3a_topX).^2)+((sourceY - intersect3a_topY)^2)+((sourceZ - intersect3a_topZ)^2));
A3b_top = sqrt(((intersect3a_topX - reflect3b_topX).^2)+((intersect3a_topY - reflect3b_topY)^2)+((intersect3a_topZ - reflect3b_topZ)^2));
B3b_top = sqrt(((intersect3b_topX - reflect3b_topX).^2)+((intersect3b_topY - reflect3b_topY)^2)+((intersect3b_topZ - reflect3b_topZ)^2));
C3b_top = sqrt(((intersect3a_topX - intersect3b_topX).^2)+((intersect3a_topY - intersect3b_topY)^2)+((intersect3a_topZ - intersect3b_topZ)^2));
A3c_top = sqrt(((reflect3c_topX - intersect3b_topX).^2)+((reflect3c_topY - intersect3b_topY)^2)+((reflect3c_topZ - intersect3b_topZ)^2));
B3c_top = sqrt(((reflect3c_topX - receiverX).^2)+((reflect3c_topY - receiverY)^2)+((reflect3c_topZ - receiverZ)^2));
C3c_top = sqrt(((receiverX - intersect3b_topX).^2)+((receiverY - intersect3b_topY)^2)+((receiverZ - intersect3b_topZ)^2));
theta3a_top = acosd((A3a_top^2 + B3a_top^2 - C3a_top^2) / (2 * A3a_top * B3a_top));
theta3b_top = acosd((A3b_top^2 + B3b_top^2 - C3b_top^2) / (2 * A3b_top * B3b_top));
theta3c_top = acosd((A3c_top^2 + B3c_top^2 - C3c_top^2) / (2 * A3c_top * B3c_top));
theta3a_top_i = theta3a_top ./ 2;
theta3b_top_i = theta3b_top ./ 2;
theta3c_top_i = theta3c_top ./ 2;
theta3a_top_t = asind((1/sqrt(3)) * sind(theta3a_top_i));
theta3b_top_t = asind((1/sqrt(3)) * sind(theta3b_top_i));
theta3c_top_t = asind((1/sqrt(3)) * sind(theta3c_top_i));
coefficient3a_top = (cosd(theta3a_top_i) - (sqrt(3) * cosd(theta3a_top_t))) / (cosd(theta3a_top_i) + (sqrt(3) * cosd(theta3a_top_t)));
coefficient3b_top = (cosd(theta3b_top_i) - (sqrt(3) * cosd(theta3b_top_t))) / (cosd(theta3b_top_i) + (sqrt(3) * cosd(theta3b_top_t)));
coefficient3c_top = (cosd(theta3c_top_i) - (sqrt(3) * cosd(theta3c_top_t))) / (cosd(theta3c_top_i) + (sqrt(3) * cosd(theta3c_top_t)));
E3_top_1 = ((1./distance3_top).*exp(-j * beta1 * distance3_top)) * (coefficient3a_top * coefficient3b_top * coefficient3c_top);
E3_top_2 = ((1./distance3_top).*exp(-j * beta2 * distance3_top)) * (coefficient3a_top * coefficient3b_top * coefficient3c_top);

% bottom
distance3_bottom = sqrt(((PtoQ - I3c_bottomX).^2)+((receiverY - I3c_bottomY)^2)+((receiverZ - I3c_bottomZ)^2));
A3a_bottom = sqrt(((sourceX - reflect3a_bottomX).^2)+((sourceY - reflect3a_bottomY)^2)+((sourceZ - reflect3a_bottomZ)^2));
B3a_bottom = sqrt(((reflect3a_bottomX - intersect3a_bottomX).^2)+((reflect3a_bottomY - intersect3a_bottomY)^2)+((reflect3a_bottomZ - intersect3a_bottomZ)^2));
C3a_bottom = sqrt(((sourceX - intersect3a_bottomX).^2)+((sourceY - intersect3a_bottomY)^2)+((sourceZ - intersect3a_bottomZ)^2));
A3b_bottom = sqrt(((intersect3a_bottomX - reflect3b_bottomX).^2)+((intersect3a_bottomY - reflect3b_bottomY)^2)+((intersect3a_bottomZ - reflect3b_bottomZ)^2));
B3b_bottom = sqrt(((intersect3b_bottomX - reflect3b_bottomX).^2)+((intersect3b_bottomY - reflect3b_bottomY)^2)+((intersect3b_bottomZ - reflect3b_bottomZ)^2));
C3b_bottom = sqrt(((intersect3a_bottomX - intersect3b_bottomX).^2)+((intersect3a_bottomY - intersect3b_bottomY)^2)+((intersect3a_bottomZ - intersect3b_bottomZ)^2));
A3c_bottom = sqrt(((reflect3c_bottomX - intersect3b_bottomX).^2)+((reflect3c_bottomY - intersect3b_bottomY)^2)+((reflect3c_bottomZ - intersect3b_bottomZ)^2));
B3c_bottom = sqrt(((reflect3c_bottomX - receiverX).^2)+((reflect3c_bottomY - receiverY)^2)+((reflect3c_bottomZ - receiverZ)^2));
C3c_bottom = sqrt(((receiverX - intersect3b_bottomX).^2)+((receiverY - intersect3b_bottomY)^2)+((receiverZ - intersect3b_bottomZ)^2));
theta3a_bottom = acosd((A3a_bottom^2 + B3a_bottom^2 - C3a_bottom^2) / (2 * A3a_bottom * B3a_bottom));
theta3b_bottom = acosd((A3b_bottom^2 + B3b_bottom^2 - C3b_bottom^2) / (2 * A3b_bottom * B3b_bottom));
theta3c_bottom = acosd((A3c_bottom^2 + B3c_bottom^2 - C3c_bottom^2) / (2 * A3c_bottom * B3c_bottom));
theta3a_bottom_i = theta3a_bottom ./ 2;
theta3b_bottom_i = theta3b_bottom ./ 2;
theta3c_bottom_i = theta3c_bottom ./ 2;
theta3a_bottom_t = asind((1/sqrt(3)) * sind(theta3a_bottom_i));
theta3b_bottom_t = asind((1/sqrt(3)) * sind(theta3b_bottom_i));
theta3c_bottom_t = asind((1/sqrt(3)) * sind(theta3c_bottom_i));
coefficient3a_bottom = (cosd(theta3a_bottom_i) - (sqrt(3) * cosd(theta3a_bottom_t))) / (cosd(theta3a_bottom_i) + (sqrt(3) * cosd(theta3a_bottom_t)));
coefficient3b_bottom = (cosd(theta3b_bottom_i) - (sqrt(3) * cosd(theta3b_bottom_t))) / (cosd(theta3b_bottom_i) + (sqrt(3) * cosd(theta3b_bottom_t)));
coefficient3c_bottom = (cosd(theta3c_bottom_i) - (sqrt(3) * cosd(theta3c_bottom_t))) / (cosd(theta3c_bottom_i) + (sqrt(3) * cosd(theta3c_bottom_t)));
E3_bottom_1 = ((1./distance3_bottom).*exp(-j * beta1 * distance3_bottom)) * (coefficient3a_bottom * coefficient3b_bottom * coefficient3c_bottom);
E3_bottom_2 = ((1./distance3_bottom).*exp(-j * beta2 * distance3_bottom)) * (coefficient3a_bottom * coefficient3b_bottom * coefficient3c_bottom);

% total
E3_1 = Ed_1 + E1_top_1 + E1_bottom_1 + E2_top_1 + E2_bottom_1 + E3_top_1 + E3_bottom_1;
E3_2 = Ed_2 + E1_top_2 + E1_bottom_2 + E2_top_2 + E2_bottom_2 + E3_top_2 + E3_bottom_2;
E3dB_1 = 20*log10(abs(E3_1));
E3dB_2 = 20*log10(abs(E3_2));

%% fourth order
% top
distance4_top = sqrt(((PtoQ - I4d_topX).^2)+((receiverY - I4d_topY)^2)+((receiverZ - I4d_topZ)^2));
A4a_top = sqrt(((sourceX - reflect4a_topX).^2)+((sourceY - reflect4a_topY)^2)+((sourceZ - reflect4a_topZ)^2));
B4a_top = sqrt(((reflect4a_topX - intersect4a_topX).^2)+((reflect4a_topY - intersect4a_topY)^2)+((reflect4a_topZ - intersect4a_topZ)^2));
C4a_top = sqrt(((sourceX - intersect4a_topX).^2)+((sourceY - intersect4a_topY)^2)+((sourceZ - intersect4a_topZ)^2));
A4b_top = sqrt(((intersect4a_topX - reflect4b_topX).^2)+((intersect4a_topY - reflect4b_topY)^2)+((intersect4a_topZ - reflect4b_topZ)^2));
B4b_top = sqrt(((intersect4b_topX - reflect4b_topX).^2)+((intersect4b_topY - reflect4b_topY)^2)+((intersect4b_topZ - reflect4b_topZ)^2));
C4b_top = sqrt(((intersect4a_topX - intersect4b_topX).^2)+((intersect4a_topY - intersect4b_topY)^2)+((intersect4a_topZ - intersect4b_topZ)^2));
A4c_top = sqrt(((reflect4c_topX - intersect4b_topX).^2)+((reflect4c_topY - intersect4b_topY)^2)+((reflect4c_topZ - intersect4b_topZ)^2));
B4c_top = sqrt(((reflect4c_topX - intersect4c_topX).^2)+((reflect4c_topY - intersect4c_topY)^2)+((reflect4c_topZ - intersect4c_topZ)^2));
C4c_top = sqrt(((intersect4c_topX - intersect4b_topX).^2)+((intersect4c_topY - intersect4b_topY)^2)+((intersect4c_topZ - intersect4b_topZ)^2));
A4d_top = sqrt(((reflect4d_topX - reflect4c_topX).^2)+((reflect4d_topY - reflect4c_topY)^2)+((reflect4d_topZ - reflect4c_topZ)^2));
B4d_top = sqrt(((receiverX - reflect4d_topX).^2)+((receiverY - reflect4d_topY)^2)+((receiverZ - reflect4d_topZ)^2));
C4d_top = sqrt(((receiverX - reflect4c_topX).^2)+((receiverY - reflect4c_topY)^2)+((receiverZ - reflect4c_topZ)^2));
theta4a_top = acosd((A4a_top^2 + B4a_top^2 - C4a_top^2) / (2 * A4a_top * B3a_top));
theta4b_top = acosd((A4b_top^2 + B4b_top^2 - C4b_top^2) / (2 * A4b_top * B3b_top));
theta4c_top = acosd((A4c_top^2 + B4c_top^2 - C4c_top^2) / (2 * A4c_top * B3c_top));
theta4d_top = acosd((A4d_top^2 + B4d_top^2 - C4d_top^2) / (2 * A4d_top * B4d_top));
theta4a_top_i = theta4a_top ./ 2;
theta4b_top_i = theta4b_top ./ 2;
theta4c_top_i = theta4c_top ./ 2;
theta4d_top_i = theta4d_top ./ 2;
theta4a_top_t = asind((1/sqrt(3)) * sind(theta4a_top_i));
theta4b_top_t = asind((1/sqrt(3)) * sind(theta4b_top_i));
theta4c_top_t = asind((1/sqrt(3)) * sind(theta4c_top_i));
theta4d_top_t = asind((1/sqrt(3)) * sind(theta4d_top_i));
coefficient4a_top = (cosd(theta4a_top_i) - (sqrt(3) * cosd(theta4a_top_t))) / (cosd(theta4a_top_i) + (sqrt(3) * cosd(theta4a_top_t)));
coefficient4b_top = (cosd(theta4b_top_i) - (sqrt(3) * cosd(theta4b_top_t))) / (cosd(theta4b_top_i) + (sqrt(3) * cosd(theta4b_top_t)));
coefficient4c_top = (cosd(theta4c_top_i) - (sqrt(3) * cosd(theta4c_top_t))) / (cosd(theta4c_top_i) + (sqrt(3) * cosd(theta4c_top_t)));
coefficient4d_top = (cosd(theta4d_top_i) - (sqrt(3) * cosd(theta4d_top_t))) / (cosd(theta4d_top_i) + (sqrt(3) * cosd(theta4d_top_t)));
E4_top_1 = ((1./distance4_top).*exp(-j * beta1 * distance4_top)) * (coefficient4a_top * coefficient4b_top * coefficient4c_top * coefficient4d_top);
E4_top_2 = ((1./distance4_top).*exp(-j * beta2 * distance4_top)) * (coefficient4a_top * coefficient4b_top * coefficient4c_top * coefficient4d_top);

% bottom
distance4_bottom = sqrt(((PtoQ - I4d_bottomX).^2)+((receiverY - I4d_bottomY)^2)+((receiverZ - I4d_bottomZ)^2));
A4a_bottom = sqrt(((sourceX - reflect4a_bottomX).^2)+((sourceY - reflect4a_bottomY)^2)+((sourceZ - reflect4a_bottomZ)^2));
B4a_bottom = sqrt(((reflect4a_bottomX - intersect4a_bottomX).^2)+((reflect4a_bottomY - intersect4a_bottomY)^2)+((reflect4a_bottomZ - intersect4a_bottomZ)^2));
C4a_bottom = sqrt(((sourceX - intersect4a_bottomX).^2)+((sourceY - intersect4a_bottomY)^2)+((sourceZ - intersect4a_bottomZ)^2));
A4b_bottom = sqrt(((intersect4a_bottomX - reflect4b_bottomX).^2)+((intersect4a_bottomY - reflect4b_bottomY)^2)+((intersect4a_bottomZ - reflect4b_bottomZ)^2));
B4b_bottom = sqrt(((intersect4b_bottomX - reflect4b_bottomX).^2)+((intersect4b_bottomY - reflect4b_bottomY)^2)+((intersect4b_bottomZ - reflect4b_bottomZ)^2));
C4b_bottom = sqrt(((intersect4a_bottomX - intersect4b_bottomX).^2)+((intersect4a_bottomY - intersect4b_bottomY)^2)+((intersect4a_bottomZ - intersect4b_bottomZ)^2));
A4c_bottom = sqrt(((reflect4c_bottomX - intersect4b_bottomX).^2)+((reflect4c_bottomY - intersect4b_bottomY)^2)+((reflect4c_bottomZ - intersect4b_bottomZ)^2));
B4c_bottom = sqrt(((reflect4c_bottomX - intersect4c_bottomX).^2)+((reflect4c_bottomY - intersect4c_bottomY)^2)+((reflect4c_bottomZ - intersect4c_bottomZ)^2));
C4c_bottom = sqrt(((intersect4c_bottomX - intersect4b_bottomX).^2)+((intersect4c_bottomY - intersect4b_bottomY)^2)+((intersect4c_bottomZ - intersect4b_bottomZ)^2));
A4d_bottom = sqrt(((reflect4d_bottomX - reflect4c_bottomX).^2)+((reflect4d_bottomY - reflect4c_bottomY)^2)+((reflect4d_bottomZ - reflect4c_bottomZ)^2));
B4d_bottom = sqrt(((receiverX - reflect4d_bottomX).^2)+((receiverY - reflect4d_bottomY)^2)+((receiverZ - reflect4d_bottomZ)^2));
C4d_bottom = sqrt(((receiverX - reflect4c_bottomX).^2)+((receiverY - reflect4c_bottomY)^2)+((receiverZ - reflect4c_bottomZ)^2));
theta4a_bottom = acosd((A4a_bottom^2 + B4a_bottom^2 - C4a_bottom^2) / (2 * A4a_bottom * B3a_bottom));
theta4b_bottom = acosd((A4b_bottom^2 + B4b_bottom^2 - C4b_bottom^2) / (2 * A4b_bottom * B3b_bottom));
theta4c_bottom = acosd((A4c_bottom^2 + B4c_bottom^2 - C4c_bottom^2) / (2 * A4c_bottom * B3c_bottom));
theta4d_bottom = acosd((A4d_bottom^2 + B4d_bottom^2 - C4d_bottom^2) / (2 * A4d_bottom * B4d_bottom));
theta4a_bottom_i = theta4a_bottom ./ 2;
theta4b_bottom_i = theta4b_bottom ./ 2;
theta4c_bottom_i = theta4c_bottom ./ 2;
theta4d_bottom_i = theta4d_bottom ./ 2;
theta4a_bottom_t = asind((1/sqrt(3)) * sind(theta4a_bottom_i));
theta4b_bottom_t = asind((1/sqrt(3)) * sind(theta4b_bottom_i));
theta4c_bottom_t = asind((1/sqrt(3)) * sind(theta4c_bottom_i));
theta4d_bottom_t = asind((1/sqrt(3)) * sind(theta4d_bottom_i));
coefficient4a_bottom = (cosd(theta4a_bottom_i) - (sqrt(3) * cosd(theta4a_bottom_t))) / (cosd(theta4a_bottom_i) + (sqrt(3) * cosd(theta4a_bottom_t)));
coefficient4b_bottom = (cosd(theta4b_bottom_i) - (sqrt(3) * cosd(theta4b_bottom_t))) / (cosd(theta4b_bottom_i) + (sqrt(3) * cosd(theta4b_bottom_t)));
coefficient4c_bottom = (cosd(theta4c_bottom_i) - (sqrt(3) * cosd(theta4c_bottom_t))) / (cosd(theta4c_bottom_i) + (sqrt(3) * cosd(theta4c_bottom_t)));
coefficient4d_bottom = (cosd(theta4d_bottom_i) - (sqrt(3) * cosd(theta4d_bottom_t))) / (cosd(theta4d_bottom_i) + (sqrt(3) * cosd(theta4d_bottom_t)));
E4_bottom_1 = ((1./distance4_bottom).*exp(-j * beta1 * distance4_bottom)) * (coefficient4a_bottom * coefficient4b_bottom * coefficient4c_bottom * coefficient4d_bottom);
E4_bottom_2 = ((1./distance4_bottom).*exp(-j * beta2 * distance4_bottom)) * (coefficient4a_bottom * coefficient4b_bottom * coefficient4c_bottom * coefficient4d_bottom);

% total
E4_1 = Ed_1 + E1_top_1 + E1_bottom_1 + E2_top_1 + E2_bottom_1 + E3_top_1 + E3_bottom_1 + E4_top_1 + E4_bottom_1;
E4_2 = Ed_2 + E1_top_2 + E1_bottom_2 + E2_top_2 + E2_bottom_2 + E3_top_2 + E3_bottom_2 + E4_top_2 + E4_bottom_2;
E4dB_1 = 20*log10(abs(E4_1));
E4dB_2 = 20*log10(abs(E4_2));

%% graph plotting
subplot(2,1,1);
hold on;
plot(PtoQ, EddB_1, 'b');
plot(PtoQ, E1dB_1, 'r');
plot(PtoQ, E2dB_1, 'k');
plot(PtoQ, E3dB_1, 'm');
plot(PtoQ, E4dB_1, 'g');
hold off;
title('Electromagnetic Field from Point p to Point q with Operation Frequency of 2.4GHz');
legend('Direct Path', 'First-Order Reflection', 'Second-Order Reflection', 'Third-Order Reflection', 'Fourth-Order Reflection');
xlabel('Line Segment "pq"/m');
ylabel('E-field/dB');

subplot(2,1,2);
hold on;
plot(PtoQ, EddB_2, 'b');
plot(PtoQ, E1dB_2, 'r');
plot(PtoQ, E2dB_2, 'k');
plot(PtoQ, E3dB_2, 'm');
plot(PtoQ, E4dB_2, 'g');
hold off;
title('Electromagnetic Field from Point p to Point q with Operation Frequency of 5.8GHz');
legend('Direct Path', 'First-Order Reflection', 'Second-Order Reflection', 'Third-Order Reflection', 'Fourth-Order Reflection');
xlabel('Line Segment "pq"/m');
ylabel('E-field/dB');

figure;
hold on;
plot(PtoQ, EddB_1, 'b');
plot(PtoQ, E1dB_1, 'r');
plot(PtoQ, E2dB_1, 'k');
plot(PtoQ, E3dB_1, 'm');
plot(PtoQ, E4dB_1, 'g');
hold off;
title('Electromagnetic Field from Point p to Point q with Operation Frequency of 2.4GHz');
legend('Direct Path', 'First-Order Reflection', 'Second-Order Reflection', 'Third-Order Reflection', 'Fourth-Order Reflection');
xlabel('Line Segment "pq"/m');
ylabel('E-field/dB');

figure;
hold on;
plot(PtoQ, EddB_2, 'b');
plot(PtoQ, E1dB_2, 'r');
plot(PtoQ, E2dB_2, 'k');
plot(PtoQ, E3dB_2, 'm');
plot(PtoQ, E4dB_2, 'g');
hold off;
title('Electromagnetic Field from Point p to Point q with Operation Frequency of 5.8GHz');
legend('Direct Path', 'First-Order Reflection', 'Second-Order Reflection', 'Third-Order Reflection', 'Fourth-Order Reflection');
xlabel('Line Segment "pq"/m');
ylabel('E-field/dB');
