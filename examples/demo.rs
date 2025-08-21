use anonyverif::*;
use std::time::Instant;

fn main() {
    // ============= Setup ================
    let mut client = Client::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER);
    let galois_keys = client.build_galois_keys();
    let keyswitch_keys = client.build_functional_keyswitch_key();
    let mut server = Server::new(DEFAULT_GLWE_PARAMTER, DEFAULT_PIR_PARAMETER, galois_keys, keyswitch_keys);

    // ============= Enroll ==============
    let db_size = 1 << 10;
    for i in 0..db_size {
        let mut features = vec![0.0; 512];
        features[i as usize % 512] = 1.0;

        client.ppclient.enroll_ggsw_masks(i as u128);
        let template_body = client.ppclient.encrypt_new_template(i, &features, 512.0);
        server.ppserver.enroll(i, template_body);
    }

    // ============== Offline Preprocess ===============
    let client_db = client.ppclient.collect_database_into_matrix().unwrap();
    let server_db = server.ppserver.collect_database_into_matrix().unwrap();

    // ============== New ct ==============
    let mut query_feature = vec![0.0f32; 512];
    query_feature[0] = 0.25;

    let glwe_ct = client.ppclient.encrypt_glwe(&query_feature, 512.0);
    let glwe_ct_cont = glwe_ct.into_container();
    let dim = glwe_ct_cont.len();
    let query_vec = nalgebra::DVector::from_iterator(dim, glwe_ct_cont.into_iter().map(|x| x >> 32));

    // ============= build query ===============
    let now = Instant::now();
    let query_idx = 0;
    let share_client = client_db.row(query_idx) * &query_vec;
    let share_client = share_client[0];
    let query_ggsw = client.construct_ggsw_monomial((share_client >> 56) as usize);
    let elapsed = now.elapsed().as_micros();
    println!("Build Query Time: {} ms, share client: {}", elapsed as f32 / 1000.0, share_client >> 56);

    // ============= build PIR query ===============
    let now = Instant::now();
    let query_ct = client.construct_query(query_idx);
    let elapsed = now.elapsed().as_micros();
    println!("Client Build Query Time: {} ms", elapsed as f32 / 1000.0);

    // ============= server PIR ================
    let now = Instant::now();
    let bs = server_db * query_vec;
    println!("similarity {}", (share_client + bs[query_idx]) >> 56);
    let mut bs = bs.into_iter().map(|&x| x as usize).collect::<Vec<_>>();
    for _ in bs.len()..(1<<16) {
        bs.push(0);
    }
    for bi in bs.iter_mut() {
        *bi >>= 56;
    }

    let pir_res = server.retrieve(&query_ct, &bs);
    let elapsed = now.elapsed().as_micros();
    println!("PIR time: {} ms", elapsed as f32 / 1000.0);

    // ============ Ext Prod ==================
    let pir_dec = client.decrypt(&pir_res);
    let now = Instant::now();
    let response_ct = server.external_product(query_ggsw, pir_res);
    let elapsed = now.elapsed().as_micros();
    println!("Extprod time: {} ms", elapsed as f32 / 1000.0);

    // =========== Dec ================
    let pt = client.decrypt(&response_ct);
    println!("res: {:?},\npir_dec: {:?}", pt, pir_dec);
}