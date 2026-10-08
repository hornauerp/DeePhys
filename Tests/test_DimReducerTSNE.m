classdef test_DimReducerTSNE < matlab.unittest.TestCase
% TEST_DIMREDUCERTSNE  t-SNE path of DimReducer (seeded, requested dimensions).
%
% Run: runtests('test_DimReducerTSNE')

    methods (Test)

        function testSeededAndDimensions(tc)
            rng(1);
            X = array2table([randn(60, 4); randn(60, 4) + 4]);
            opts = struct('Method', 'tSNE', 'NDims', 3, 'Seed', 7);
            r1 = DimReducer.reduce(X, opts);
            r2 = DimReducer.reduce(X, opts);
            tc.verifyEqual(size(r1.Reduction), [120 3]);
            tc.verifyEqual(r1.Reduction, r2.Reduction);   % same seed, same embedding
        end

    end
end
