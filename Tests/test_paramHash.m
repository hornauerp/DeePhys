classdef test_paramHash < matlab.unittest.TestCase
% TEST_PARAMHASH  paramHash / md5hex: identical IDs with and without a JVM.
%
% Run: runtests('test_paramHash')

    methods (Test)

        function testRFC1321Vectors(tc)
            tc.verifyEqual(md5hex(uint8('')),    'd41d8cd98f00b204e9800998ecf8427e');
            tc.verifyEqual(md5hex(uint8('abc')), '900150983cd24fb0d6963f7d28e17f72');
            tc.verifyEqual(md5hex(uint8('message digest')), 'f96b697d7cb7938d525a2f31aaf161d0');
            tc.verifyEqual(md5hex(uint8(['12345678901234567890123456789012345678901234567890' ...
                '123456789012345678901234567890'])), '57edf4a22be3c955ac49da2e2107b67a');
        end

        function testMatchesJavaAcrossPaddingBoundaries(tc)
            tc.assumeTrue(usejava('jvm'), 'Needs a JVM for the reference digest');
            rng(3);
            for n = [1 55 56 57 63 64 65 119 120 1000 4099]
                bytes = uint8(randi([0 255], n, 1));
                md = java.security.MessageDigest.getInstance('MD5');
                md.update(bytes);
                ref = sprintf('%02x', typecast(md.digest(), 'uint8'));
                tc.verifyEqual(md5hex(bytes), ref, sprintf('length %d', n));
            end
        end

        function testParamHashIsDeterministic(tc)
            s = struct('ChipID', "CHIP0001_12", 'PlatingDate', 240101, 'TemplateID', 17);
            tc.verifyEqual(paramHash(s), paramHash(s));
            tc.verifyNotEqual(paramHash(s), paramHash(setfield(s, 'TemplateID', 18))); %#ok<SFLD>
            tc.verifyEqual(strlength(string(paramHash(s))), 32);
        end

    end
end
