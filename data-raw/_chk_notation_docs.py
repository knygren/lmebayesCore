import io, re
for p in ['inst/notation.md',
          'inst/HIERARCHICAL_LINEAR_MODEL_NOTATION.md',
          'inst/HIERARCHICAL_GENERALIZED_LINEAR_MODEL_NOTATION.md']:
    s = io.open(p, encoding='utf-8').read()
    body = re.sub(r'\$\$.*?\$\$', '', s, flags=re.S)
    print('%-52s  $$:%3d  $:%4d  braces:%3d  headings:%2d  lines:%4d'
          % (p.split('/')[-1], s.count('$$'), body.count('$'),
             s.count('{') - s.count('}'),
             len(re.findall(r'^## ', s, flags=re.M)), s.count('\n')))
    for i, line in enumerate(s.split('\n'), 1):
        if line.startswith('|') and line.count('$') % 2:
            print('   odd $ count, line %d: %s' % (i, line[:70]))
