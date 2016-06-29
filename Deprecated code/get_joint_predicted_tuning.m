function [something] = get_joint_predicted_tuning(bdf,workspace,joint_weights)
% GET_JOINT_PREDICTED_TUNING Get predicted tuning, assuming neurons draw
% directly from joint inputs. Assumes GLM for generative model and tuning.
%   Inputs:
%       BDF - postprocessed bdf, including OpenSim kinematics as joint_pos
%       and muscle_pos
%       WORKSPACE - pair of numbers indicating corners of workspace. Order
%       is [far-left close-right] or [top-left bottom-right]
%       JOINT_WEIGHTS - weights on 7 joints. Includes first column of
%       baseline values.