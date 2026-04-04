function [loaded_ok, DCM] = load_dcm_and_psd(full_path, filename)
    loaded_ok = false;
    DCM = struct();

    loaded_data = load(full_path, 'DCM');
    if ~isfield(loaded_data, 'DCM')
        fprintf('      WARNING: File does not contain DCM field: %s\n', filename);
        return;
    end

    DCM = loaded_data.DCM;
    if ~(isfield(DCM, 'xY') && isfield(DCM.xY, 'y') && ~isempty(DCM.xY.y))
        fprintf('      WARNING: DCM lacks proper data structure in %s\n', filename);
        return;
    end

    loaded_ok = true;
end
