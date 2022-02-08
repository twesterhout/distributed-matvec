module Distribute {
  use BlockDist;
  // use BlockCycDist;
  // use CyclicDist;
  use states only hash64_01;
  use basis only BasisStates;

  proc distributionMask(const ref states : [?D] uint(64)) {
    const n = numLocales;
    var mask : [D] int;
    forall (x, i) in zip(states, mask) {
      assert(x.locale == here);
      assert(i.locale == here);
      i = (hash64_01(x) % n:uint):int;
    }
    return mask;
  }

  proc distributionMask(const ref states : BasisStates, chunkSize : int) {
    const totalCount = + reduce states.counts;
    const D : domain(1) dmapped Block(LocaleSpace) = {0 ..# totalCount};
    var mask : [D] int;
    forall (offset, indices) in kMergeIndicesChunked(states.representatives,
                                                     states.counts, chunkSize) {
      mask[offset ..# indices.size] = [i in indices.domain] indices[i][0];
    }
    // writeln(D);
    // for i in LocaleSpace {
    //   writeln(mask[mask.localSubdomain(Locales[i])]);
    // }
    return mask;
  }

  private proc distributionCounts(const ref mask : [] int) {
    var histogram : [LocaleSpace] int;
    forall i in mask with (+ reduce histogram) {
      histogram[i] += 1;
    }
    return histogram;
  }

  private proc maxInnerCount(const ref mask : [] int)
      where mask.hasSingleLocalSubdomain() {
    var globalCounts : [LocaleSpace] [LocaleSpace] int;
    coforall loc in Locales do on loc {
      const localIndices = mask.localSubdomain();
      const localCounts = distributionCounts(mask[localIndices]);
      globalCounts[loc.id] = localCounts;
    }
    return max reduce ([i in LocaleSpace] max reduce globalCounts[i]);
  }

  private inline proc innerDomain(const ref xs : [?D] ?eltType, n : int) : domain(1)
      where D.rank == 1 {
    return {0 ..# n};
  }
  private inline proc innerDomain(const ref xs : [?D] ?eltType, n : int) : domain(2)
      where D.rank == 2 {
    return {0 ..# D.dim(0).size, 0 ..# n};
  }

  proc distributeArray(const ref array : [?OuterDomain] ?eltType,
                       const ref mask  : [?OuterDomain2] int)
      where isSubtype(OuterDomain.dist.type, Block) {
    assert(array.size > 0);
    param rank = array.domain.rank;
    const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);

    const TempInnerDomain = innerDomain(array, maxInnerCount(mask));
    var globalTemp : [OnePerLocale] [LocaleSpace] [TempInnerDomain] eltType;
    var globalTempCounts : [OnePerLocale] [LocaleSpace] int;
    coforall loc in Locales do on loc {
      ref temp = globalTemp[loc.id];
      ref offsets = globalTempCounts[loc.id];
      for _i in mask.localSubdomain() do {
        const targetLocale = mask[_i];
        if rank == 1 {
          temp[targetLocale][offsets[targetLocale]] = array[_i];
        }
        else if rank == 2 {
          for k in TempInnerDomain.dim(0) {
            temp[targetLocale][k, offsets[targetLocale]] = array[k, _i];
          }
        }
        else { compilerError("Oops"); }
        offsets[targetLocale] += 1;
      }
    }

    const globalFinalCounts : [OnePerLocale] int =
      [j in LocaleSpace] (+ reduce [i in OnePerLocale] globalTempCounts[i][j]);
    const maxCountPerLocale = max reduce globalFinalCounts;
    const FinalInnerDomain = innerDomain(array, maxCountPerLocale);
    var globalFinal : [OnePerLocale] [FinalInnerDomain] eltType;
    coforall loc in Locales do on loc {
      var offset : int = 0;
      for i in LocaleSpace {
        const n = globalTempCounts[i][loc.id];
        if rank == 1 {
          globalFinal[loc.id][offset ..# n] = globalTemp[i][loc.id][0 ..# n];
        }
        else if rank == 2 {
          globalFinal[loc.id][.., offset ..# n] = globalTemp[i][loc.id][.., 0 ..# n];
        }
        offset += n;
      }
    }

    return (globalFinal, globalFinalCounts);
  }

  inline proc distributeStates(const ref states) {
    return distributeStates(states, distributionMask(states));
  }
  proc distributeStates(const ref states : [?D] uint(64),
                        const ref mask : [?D2] int)
      where D.rank == 1 {
    var (chunks, counts) = distributeArray(states, mask); 
    return new BasisStates(chunks[0].size, chunks, counts);
  }

  proc distributeVectors(const ref vectors : [?D] ?eltType,
                         const ref mask : [?D2] int)
      where D.rank == 2 {
    var (chunks, counts) = distributeArray(vectors, mask); 
    return chunks;
  }
} // module Distribute
