function write_summary_csv(csv_path, header, rows)
    fid = fopen(csv_path, 'w');
    if fid < 0
        error('Could not open CSV for writing: %s', csv_path);
    end
    cleanup_obj = onCleanup(@() fclose(fid));

    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', header{:});
    for i = 1:size(rows, 1)
        stage = rows{i, 1};
        receptor = rows{i, 2};
        n_adj = rows{i, 3};
        adj = rows{i, 4};
        mse = rows{i, 5};
        es = rows{i, 6};
        pval = rows{i, 7};
        qval = rows{i, 8};
        baseline_source = rows{i, 9};
        best_file = rows{i, 10};
        simdig = rows{i, 11};
        fprintf(fid, '%s,%s,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%s,%s,%.12g\n', ...
            stage, receptor, n_adj, adj, mse, es, pval, qval, baseline_source, best_file, simdig);
    end
end
