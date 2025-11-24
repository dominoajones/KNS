step = []; % Initialize empty array for storing steps
tol = 0.1; % Define a small tolerance for floating-point comparison

% Loop through steps in TransientSolution
for i = 1:(md.results.TransientSolution(end).step / md.settings.output_frequency)
   % Check if the time at the current step matches any value in t
   for p = 1:length(t)
      if abs(md.results.TransientSolution(i).time - t(p)) < tol
         step = [step, i]; % Append the index to step
         break; % Exit inner loop since a match is found
      end
   end
end

