function filtered_data = filter_patient_data(data)
    filtered_data = {};
    for i = 1:length(data)
        p = dcm_wake_parse(data{i});
        if p.ok && p.is_patient
            filtered_data{end+1} = data{i};
        end
    end
    filtered_data = filtered_data';
end
