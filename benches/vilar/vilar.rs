#![allow(mixed_script_confusables, unused_variables)]

use rebop::define_system;

define_system! {
    αA αpA αR αpR βA βR δMA δMR δA δR γA γR γC θA θR;
    Vilar { Da, Dr, Dpa, Dpr, Ma, Mr, A, R, C }
    r_activation_a      : Da + A    => Dpa      @ γA
    r_activation_r      : Dr + A    => Dpr      @ γR
    r_deactivation_a    : Dpa       => Da + A   @ θA
    r_deactivation_r    : Dpr       => Dr + A   @ θR
    r_transcription_a   : Da        => Da + Ma  @ αA
    r_transcription_r   : Dr        => Dr + Mr  @ αR
    r_transcription_p_a : Dpa       => Dpa + Ma @ αpA
    r_transcription_p_r : Dpr       => Dpr + Mr @ αpR
    r_translation_a     : Ma        => Ma + A   @ βA
    r_translation_r     : Mr        => Mr + R   @ βR
    r_complexation      : A + R     => C        @ γC
    r_decomplexation    : C         => R        @ δA
    r_decay_mRNA_a      : Ma        =>          @ δMA
    r_decay_mRNA_r      : Mr        =>          @ δMR
    r_decay_prot_a      : A         =>          @ δA
    r_decay_prot_r      : R         =>          @ δR
}

fn main() {
    for _ in 0..100 {
        let mut vilar = Vilar::new();
        vilar.αA = 50.;
        vilar.αpA = 500.;
        vilar.αR = 0.01;
        vilar.αpR = 50.;
        vilar.βA = 50.;
        vilar.βR = 5.;
        vilar.δMA = 10.;
        vilar.δMR = 0.5;
        vilar.δA = 1.;
        vilar.δR = 0.2;
        vilar.γA = 1.;
        vilar.γR = 1.;
        vilar.γC = 2.;
        vilar.θA = 50.;
        vilar.θR = 100.;
        vilar.Da = 1;
        vilar.Dr = 1;
        for i in 1..=200 {
            vilar.advance_until(i as f64);
        }
        println!("{}, {}, {}", vilar.A, vilar.R, vilar.C);
    }
}
