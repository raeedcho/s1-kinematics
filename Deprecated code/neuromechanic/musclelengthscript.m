%%
nm = xml2matlab('monkey.nmco');                  % Convert xml to a matlab cell structure
bod = nm.NeuromechanicFile{1}.Bodies{1}.RigidBody;  % Get <RigidBody> elements
mass = 0;                                           % Initialize mass to zero
for ii = 1:length(bod)
    mass = mass+bod{ii}.Mass{1}.Value;              % Add total mass from each body
end;
%%
dyn = nm.NeuromechanicFile{1}.Dynamics;             % Get <Dynamic> elements
n = length(dyn);                                    % Get the number of data points
t = zeros(n,1);                                     % Initialize time
musc_len = zeros(n,39);
for ii = 1:n
    t(ii,1) = dyn{ii}.Time{1}.Value;                % Store time
    musc_len(ii,:) = dyn{ii}.MuscleLength{1}.Value;
end;
%%
plot(t,musc_len)
xlabel('Time (sec)')
ylabel('Muscle length (m)')