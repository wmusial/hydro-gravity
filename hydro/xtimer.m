classdef xtimer

  properties (GetAccess='public', SetAccess='private')
    tbeg = containers.Map;
    ttot = containers.Map;
    tstart = 0;
  end

  methods
    function this=xtimer()
      this.tstart = tic;
      this.tbeg = containers.Map;
      this.ttot = containers.Map;
    end;

    function tic(this, tname)
      this.ensureTimer(tname);

      this.tbeg(tname) = tic;
    end;

    function toc(this, tname)
      if ~this.tbeg.isKey(tname)
        return
      end

      if this.tbeg(tname) == 0
        return
      end

      this.ensureTimer(tname);
      this.ttot(tname) = this.ttot(tname) + toc(this.tbeg(tname));
    end

    function r = ensureTimer(this, tname)
      r = false;
      if this.ttot.isKey(tname)
        r = true;
      else
        this.ttot(tname) = 0;
        this.tbeg(tname) = 0;
      end
    end

    function stats(this)
      fprintf('\n***** timer statistics *****\n');
      vals = this.ttot.values;
      keys = this.ttot.keys;

      tot = toc(this.tstart);
      ttot = sum([vals{:}]);

      fprintf('total time:            %.1fs\n', tot);
      fprintf('total registered time: %.1fs\n\n', ttot);

      tname_length = cellfun(@(x)(length(x)), keys);
      tname_length_max = max(tname_length);

      ttot_length = cellfun(@(x)(length(num2str(round(x)))), vals);
      ttot_length_max = max(ttot_length);

      for i=1:length(keys)
        fprintf(['%' num2str(tname_length_max) 's : %' num2str(ttot_length_max+2) '.1fs  %5.1f%%  %5.1f%%\n'], keys{i}, vals{i}, vals{i} ./ ttot * 100, vals{i} ./ tot * 100);
      end
    end
  end
end
