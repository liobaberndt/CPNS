function parsave(filename, variable)
    try
        sim_psd = variable;
        save(filename, 'sim_psd', '-mat');
        fprintf('    Successfully saved file: %s\n', filename);
    catch ME
        fprintf('    ERROR saving file %s: %s\n', filename, ME.message);
    end
end
