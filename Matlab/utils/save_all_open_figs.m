function save_all_open_figs(outdir, base)
    if ~exist(outdir,"dir"), mkdir(outdir); end
    figs = findall(groot,'Type','figure');
    for k = 1:numel(figs)
        nm = string(figs(k).Name);
        if strlength(nm)==0, nm = base + "_" + sprintf("%02d",k); end
        exportgraphics(figs(k), fullfile(outdir, nm + ".png"), "Resolution",300);
        exportgraphics(figs(k), fullfile(outdir, nm + ".pdf"), "ContentType","vector");
        savefig(figs(k), fullfile(outdir, nm + ".fig"));
    end
end