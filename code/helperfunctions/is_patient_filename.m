function is_patient = is_patient_filename(filename)
    p = dcm_wake_parse(filename);
    is_patient = p.ok && p.is_patient;
end
