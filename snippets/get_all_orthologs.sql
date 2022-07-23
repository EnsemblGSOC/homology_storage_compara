SET @seqid = (SELECT canonical_member_id FROM gene_member WHERE stable_id = 'ENSG00000182463')
SELECT * FROM homology_member
WHERE homology_id IN (
    SELECT hm.homology_id
    FROM homology_member hm LEFT JOIN homology h ON h.homology_id = hm.homology_id
    WHERE hm.seq_member_id = @seqid AND description LIKE 'ortholog_one2one'
) AND seq_member_id <> @seqid
UNION
SELECT * FROM homology_member
WHERE homology_id IN (
    SELECT hm.homology_id
    FROM homology_member hm LEFT JOIN homology h ON h.homology_id = hm.homology_id
    WHERE hm.seq_member_id = @seqid AND description LIKE 'ortholog_many2many'
) AND seq_member_id <> @seqid;