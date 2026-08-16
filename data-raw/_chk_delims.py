import re, io
p = 'inst/RESTRICTED_GIBBS_MINORIZATION_TV.md'
s = io.open(p, encoding='utf-8').read()
for a, b in [(r'\[', r'\]'), (r'\(', r'\)')]:
    na = len(re.findall(re.escape(a), s))
    nb = len(re.findall(re.escape(b), s))
    print('%-4s %5d   %-4s %5d   %s' % (a, na, b, nb, 'OK' if na == nb else 'MISMATCH'))
print('boxed  %5d' % s.count(r'\boxed'))
print('braces %5d (want 0)' % (s.count('{') - s.count('}')))
print('lines  %5d' % s.count('\n'))
for tag in ['2A.3-G', 'S_\\flat', 'Lemma 2A.3-F', 'Kuratowski']:
    print('%-14s %4d' % (tag, s.count(tag)))
