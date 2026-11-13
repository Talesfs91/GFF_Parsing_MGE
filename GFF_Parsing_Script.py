import os
import sys

def parse_gff_native(gff_file):
    """
    Parses a GFF3 file natively without third-party dependencies.
    """
    targets, mges = [], []
    mge_keywords = ['transposase', 'integrase', 'insertion sequence', 'is', 'integron', 'ice', 'relaxase']
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            seq_id, _, feat_type, start, end, _, strand, _, attributes = parts
            
            # Extract attributes (ID, Name, Product, Note, Gene)
            attr_dict = {}
            for item in attributes.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k.lower()] = v
            
            name = attr_dict.get('gene', attr_dict.get('name', attr_dict.get('id', 'Unknown')))
            product = attr_dict.get('product', '').lower()
            note = attr_dict.get('note', '').lower()
            
            feature_dict = {
                'seq_id': seq_id,
                'name': name,
                'product': attr_dict.get('product', 'N/A'),
                'start': int(start),
                'end': int(end),
                'strand': strand
            }
            
            # Classify as Mobile Genetic Element or regular Gene/CDS
            if any(kw in product or kw in note or kw in feat_type.lower() for kw in mge_keywords):
                mges.append(feature_dict)
            elif feat_type in ['gene', 'CDS']:
                targets.append(feature_dict)
                
    return targets, mges

def analyze_distances(targets, mges, max_distance=10000):
    """
    Calculates distances from target genes to MGEs and checks for flanking setups.
    """
    results = []
    
    targets = sorted(targets, key=lambda x: (x['seq_id'], x['start']))
    mges = sorted(mges, key=lambda x: (x['seq_id'], x['start']))
    
    for target in targets:
        seq_id = target['seq_id']
        t_start, t_end = target['start'], target['end']
        
        local_mges = [m for m in mges if m['seq_id'] == seq_id]
        
        upstream_mges = [m for m in local_mges if m['end'] <= t_start]
        downstream_mges = [m for m in local_mges if m['start'] >= t_end]
        
        closest_upstream = upstream_mges[-1] if upstream_mges else None
        closest_downstream = downstream_mges[0] if downstream_mges else None
        
        dist_upstream = (t_start - closest_upstream['end']) if closest_upstream else float('inf')
        dist_downstream = (closest_downstream['start'] - t_end) if closest_downstream else float('inf')
        
        is_flanked = (dist_upstream <= max_distance) and (dist_downstream <= max_distance)
        min_dist = min(dist_upstream, dist_downstream)
        
        if min_dist <= max_distance:
            status = "HIGH POTENTIAL (FLANKED)" if is_flanked else "POTENTIALLY MOBILE"
            
            results.append({
                'gene': target['name'],
                'product': target['product'],
                'contig': seq_id,
                'min_distance_bp': min_dist,
                'upstream_mge': closest_upstream['name'] if closest_upstream else 'None',
                'upstream_distance_bp': dist_upstream if closest_upstream else None,
                'downstream_mge': closest_downstream['name'] if closest_downstream else 'None',
                'downstream_distance_bp': dist_downstream if closest_downstream else None,
                'status': status
            })
            
    return results

if __name__ == "__main__":
    gff_filename = sys.argv[1]
    MAX_WINDOW = 10000 
    
    if not os.path.exists(gff_filename):
        print(f"Error: File '{gff_filename}' not found.")
    else:
        genes, mges = parse_gff_native(gff_filename)
        print(f"Found {len(genes)} target genes and {len(mges)} Mobile Genetic Elements.")
        
        print(f"\nAnalyzing proximity window ({MAX_WINDOW/1000} kb)...")
        hgt_candidates = analyze_distances(genes, mges, max_distance=MAX_WINDOW)
        
        print("\n--- HGT TRANSFER POTENTIAL REPORT ---")
        print("Gene\tMin Dist (bp)\tUp MGE\tUp Dist (bp)\tDown MGE\tDown Dist (bp)\tStatus")
        print("-" * 120)
        for cand in hgt_candidates:
            up_dist = cand['upstream_distance_bp'] if cand['upstream_distance_bp'] is not None else 'None'
            down_dist = cand['downstream_distance_bp'] if cand['downstream_distance_bp'] is not None else 'None'
            print(f"{cand['gene']}\t{cand['min_distance_bp']}\t{cand['upstream_mge']}\t{up_dist}\t{cand['downstream_mge']}\t{down_dist}\t{cand['status']}")