#모듈 불러오기
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter, NullFormatter # <-- 이 줄이 빠져서 발생하는 에러입니다.
import streamlit as st
import pandas as pd
import io
import os
import ssl
import subprocess
import requests
import matplotlib.pyplot as plt
import numpy as np
import json
from Bio.Seq import Seq  # siRNA 역상보 서열 계산용
from Bio.SeqUtils import MeltingTemp as mt
import free_energy
import general_helpers

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML

from stmol import showmol
import py3Dmol


#보안 및 기본 설정
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "csw1567@naver.com" 

if 'setup_done' not in st.session_state:
    st.set_page_config(
        page_title="RNAi design prototype(KNU)",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    st.session_state['setup_done'] = True




#모든 텝에 적용될 CSS 디자인 넣기
st.markdown("""
    <style>
    /* 1. 탭 전체 글자 크기 및 스타일 */
    .stTabs [data-baseweb="tab-list"] button p {
        font-size: 18px;
        font-weight: bold;
    }

    /* 2. 탭 버튼 사이에 세로선 넣기 */
    .stTabs [data-baseweb="tab-list"] button {
        border-right: 1px solid #d3d3d3; /* 오른쪽에 회색 세로선 추가 */
        border-radius: 0px; /* 각진 모양으로 변경하여 선이 깔끔하게 보이게 함 */
        padding-left: 20px;
        padding-right: 20px;
    }

    /* 3. 마지막 탭은 세로선 없애기 */
    .stTabs [data-baseweb="tab-list"] button:last-child {
        border-right: none;
    }

    /* 4. 선택된 탭 배경색 및 강조 */
    .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
        background-color: #e6f3ff;
        color: #1a2a6c;
    }
    </style>
    """, unsafe_allow_html=True)

# 4. 사이드바 구성 
with st.sidebar:
    st.title("KNU Bio-Lab")
    st.caption("Insect Physiology Laboratory")
    st.markdown("---")

    # 문의처 섹션 (HTML/CSS 사용)
    st.markdown("""
        <div style="background-color:#f0f2f6; padding:15px; border-radius:10px; border-left:5px solid #1a2a6c; margin-bottom:20px;">
            <div style="color:#1a2a6c; font-weight:bold; margin-bottom:5px;">📧 Inquiry & Support</div>
            <p style="margin:0; font-size:0.95rem;"><b>곤충생리학연구실</b></p>
            <p style="margin:0; font-size:0.9rem;">최성환 연구원에게 찾아오세요.</p>
        </div>
    """, unsafe_allow_html=True)
    
    # 유튜브 영상 섹션 
    st.markdown("### 📺 RNAi Learning Resource")
    st.video("https://www.youtube.com/watch?v=zVe3Nom2SpI")
    st.caption("Mechanism of RNA Interference")
    st.markdown("---")

    
    
    # 버전 정보 (맨 밑)
    st.markdown("""
        <div style="font-size:0.8rem; color:grey; text-align:center;">
            System Version: 1.1.0-beta<br>
            © 2026 KNU Insect Physiology Lab
        </div>
    """, unsafe_allow_html=True)


#메인 타이틀 디자인
st.title("Bursaphelenchus xylophilus genetic analysis tool")
tab1, tab2, tab3, tab4, tab5, tab6,  = st.tabs([
    "타겟 유전자 정보 조회",
    "RNAi디자인 및 off-target",
    "가상 gel imagine 생성",
    "농도 계산기",
    "Gene visualization",
    "Gene조정 및 파일형식 변환"
])

import re
import urllib.parse
import subprocess
import tempfile
import re





# 2. 함수 정의 (호출보다 위에 있어야 함)
def normalize_id(full_id):
    if full_id is None or str(full_id) == 'nan':
        return ""
    # 데이터를 강제로 문자로 바꾸고 분리
    return str(full_id).split('.')[0].strip()



if os.path.exists(result_csv) and os.path.getsize(result_csv) > 0:
    # 3. CSV 읽을 때 컬럼명 명시
    df = pd.read_csv(result_csv, names=["Query", "Locus ID", "Identity(%)", "Length", "Mismatch", "Gaps", "Q_Start", "Q_End", "S_Start", "S_End", "E-value", "BitScore"])
    
    # 4. 안전하게 적용 (astype(str) 추가)
    if 'Locus ID' in df.columns:
        df['Normalized ID'] = df['Locus ID'].astype(str).apply(normalize_id)
        
        # 이름 매핑 (get_descriptions()가 사전 데이터를 반환한다고 가정)
        names_dict = get_descriptions()
        df['Target Function'] = df['Normalized ID'].map(lambda x: names_dict.get(x, "No Map Found"))
        
        st.dataframe(df[["Target Function", "Normalized ID", "Identity(%)"]])
with tab1:
    st.header("primer를 통한 target찾기")

    # [1. 이름표 로드: 이제 ID가 100% 일치하므로 정규식 필요 없음]
    @st.cache_data
    def get_descriptions():
        desc_dict = {}
        current_dir = os.path.dirname(os.path.abspath(__file__))
        target_path = os.path.join(current_dir, "pwn_pro_named.fa")
        
        if os.path.exists(target_path):
            from Bio import SeqIO
            for record in SeqIO.parse(target_path, "fasta"):
                full_desc = record.description
                # ">ID 이름" 형태에서 이름만 추출
                parts = full_desc.split(" ", 1)
                name = parts[1] if len(parts) > 1 else "Hypothetical Protein"
                desc_dict[str(record.id)] = name
        return desc_dict

    # [2. UI 섹션]
    query_seq = st.text_area("프라이머 서열 입력", height=100, 
                             placeholder="예: GCTTCACAGCAGTGCCAAG", key="final_system")

    if st.button("타겟 유전자 분석 실행", use_container_width=True):
        if query_seq:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # 새로 만든 DB 경로
            db_path = os.path.join(current_dir, "pwn_fiXed_db", "pwn_fixed_db")
            temp_query = os.path.join(current_dir, "temp_query.fa")
            result_csv = os.path.join(current_dir, "blast_result.csv")

            with st.spinner("재선충 전용 DB 검색 중..."):
                with open(temp_query, "w") as f:
                    f.write(f">Query\n{query_seq}")

                import subprocess
                import pandas as pd
                
                # BLAST 실행
                subprocess.run(["blastn", "-query", temp_query, "-db", db_path, "-out", result_csv,
                               "-outfmt", "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                               "-task", "blastn-short", "-evalue", "1000"])

                if os.path.exists(result_csv) and os.path.getsize(result_csv) > 0:
                    df = pd.read_csv(result_csv, names=["Query", "Locus ID", "Identity(%)", "Length", "Mismatch", "Gaps", "Q_Start", "Q_End", "S_Start", "S_End", "E-value", "BitScore"])
                    
                    names = get_descriptions()
                    
                    # [매핑] 이제 1:1로 완벽하게 매칭됩니다!
                    df['Target Function'] = df['Locus ID'].apply(lambda x: names.get(str(x), "No Map Found"))
                    df['NCBI Link'] = df['Locus ID'].apply(lambda x: f"https://www.ncbi.nlm.nih.gov/search/all/?term={x}")

                    st.success(f"분석 완료! 타겟 유전자를 확인하세요.")
                    st.dataframe(
                        df[["Target Function", "Locus ID", "NCBI Link", "Identity(%)", "E-value"]],
                        column_config={"NCBI Link": st.column_config.LinkColumn("NCBI", display_text="Search ")},
                        use_container_width=True, hide_index=True
                    )
                else:
                    st.error("매칭되는 유전자가 없습니다. 서열을 확인해 주세요.")
with tab2:
    st.header("🧬 si-Fi RNAi 분석 엔진")
    st.info("CDS 파일을 업로드하여 최적의 siRNA 후보군을 탐색합니다.")

    # 1. 파일 업로드 칸 추가
    st.markdown("### 1. 시퀀스 파일 업로드")
    uploaded_file = st.file_uploader("분석할 CDS 또는 FASTA 파일을 선택하세요", type=["fasta", "fa", "cds"])
    
    # 텍스트 입력도 병행하고 싶을 경우를 위한 처리
    target_fasta = ""
    if uploaded_file is not None:
        # 업로드된 파일 내용을 읽어옴
        target_fasta = uploaded_file.read().decode("utf-8")
        st.success(f"파일 업로드 완료: {uploaded_file.name}")
    else:
        # 파일이 없을 때만 직접 입력창을 보여줌
        target_fasta = st.text_area("또는 서열을 직접 입력하세요 (FASTA 형식)", height=150)

    # 2. 파라미터 설정 (기존 코드 유지)
    st.markdown("---")
    col_mode, col_set = st.columns(2)
    with col_mode:
        mode_selection = st.radio("분석 모드", ["RNAi Design Mode", "Off-target Prediction Mode"])
        mode = 0 if mode_selection == "RNAi Design Mode" else 1
    with col_set:
        db_name = st.text_input("데이터베이스 이름", value="pwn_db")
        si_size = st.number_input("siRNA 크기 (bp)", value=21)

    # 3. 실행 로직
    if st.button("🚀 siRNA 후보군 분석 시작", use_container_width=True):
        if not target_fasta.strip():
            st.error("파일을 업로드하거나 서열을 입력해야 합니다.")
        else:
            work_dir = r"C:\Users\SAMSUNG\blast_work"
            
            with st.spinner("CDS 분석 및 후보군 추출 중..."):
                try:
                    # 입력 서열 임시 파일 저장
                    tmp_query = os.path.join(work_dir, "query_tmp.fa")
                    with open(tmp_query, "w", encoding="utf-8") as f:
                        f.write(target_fasta.strip())

                    # 파이프라인 호출
                    pipeline = SifiPipeline(
                        bowtie_db=db_name,
                        db_location=work_dir + os.sep,
                        query_sequences=tmp_query,
                        sirna_size=si_size,
                        # ... (나머지 인자들은 기존과 동일하게 전달)
                    )

                    results = pipeline.run_pipeline
                    
                    if results:
                        img, json_path, table_data, main_dict, off_dict, err = results
                        
                        if err:
                            st.error(f"오류: {err}")
                        else:
                            st.success("✅ 분석이 완료되었습니다. 아래에서 후보군을 확인하세요.")
                            # 결과 출력 로직...
                            
                except Exception as e:
                    st.error(f"에러 발생: {e}")

with tab3:
    st.header("Multi-Lane Virtual Gel Simulator")
    st.info("각 bp 값을 쉼표(,)로 구분해 입력하면 순서대로 Lane 1, 2, 3...에 배치됩니다.")

    # 1. 입력 섹션
    raw_bp_inputs = st.text_input(
        "Enter Band Sizes for each Lane (e.g., 500, 1200, 800)", 
        value="500, 1200, 800",
        key="gel_input_field"
    )
    
    # --- 변수 초기화 및 입력 처리 (에러 방지 핵심 구간) ---
    target_bp_list = []
    if raw_bp_inputs:
        try:
            # 쉼표로 나누고 공백 제거 후 정수 변환
            target_bp_list = [int(x.strip()) for x in raw_bp_inputs.split(',') if x.strip()]
        except ValueError:
            st.error("⚠️ 숫자 형식이 올바르지 않습니다. 숫자와 쉼표만 사용해 주세요 (예: 500, 1000).")
            target_bp_list = []

    # 2. 결과 시각화 섹션 (target_bp_list가 존재할 때만 실행)
    if target_bp_list:
        with st.spinner("Generating Detailed Gel..."):
            # 레인 개수 설정 (Ladder 포함 최소 5개 레인 생성)
            num_lanes = max(5, len(target_bp_list) + 1) 
            fig, ax = plt.subplots(figsize=(num_lanes * 1.5, 9)) 
            ax.set_facecolor('#361F00') # 사진과 같은 짙은 브라운 배경
            
            # Y축 로그 스케일 및 범위 설정
            ax.set_yscale('log')
            ax.set_ylim(40, 6000)
            ax.set_xlim(0.5, num_lanes + 0.5)

            # --- 눈금자 설정 (모든 밴드 위치에 숫자 표시) ---
            ladder_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000]
            
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
            
            ax.set_yticks(ladder_sizes) 
            ax.set_yticklabels([f"{x}" for x in ladder_sizes], color='#FFA54F', fontsize=8)
            
            # X축 라벨 설정
            lane_labels = ["L"] + [f"Lane {i+1}" for i in range(num_lanes-1)]
            ax.set_xticks(range(1, num_lanes + 1))
            ax.set_xticklabels(lane_labels, color='#FFA54F', fontweight='bold')
            
            # 눈금 및 테두리 스타일
            ax.tick_params(axis='y', colors='#FFA54F', length=4, labelsize=8)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.grid(False)

            # --- 메커니즘 1: 상단 Well (시료 주입구) ---
            for x in range(1, num_lanes + 1):
                well = Rectangle((x - 0.35, 4800), 0.7, 800, color='#D2691E', alpha=0.7)
                ax.add_patch(well)

            # --- 메커니즘 2: Lane L (Ladder) 밴드 그리기 ---
            for l_size in ladder_sizes:
                # 주요 지점(500, 1000, 3000)은 더 두껍게
                lw = 2.5 if l_size in [500, 1000, 3000] else 1.2
                ax.hlines(y=l_size, xmin=0.6, xmax=1.4, colors='#FF7F24', linewidth=lw, alpha=0.7)

            # --- 메커니즘 3: 각 Lane에 입력된 밴드 배정 ---
            for i, bp in enumerate(target_bp_list):
                lane_idx = i + 2 # Ladder가 1번이므로 2번 레인부터 시작
                if lane_idx <= num_lanes:
                    # 메인 밴드 (강한 주황 형광)
                    ax.hlines(y=bp, xmin=lane_idx - 0.35, xmax=lane_idx + 0.35, 
                              colors='#FF7F24', linewidth=6, alpha=0.95)
                    # 번짐 효과 (Glow)
                    ax.hlines(y=bp, xmin=lane_idx - 0.4, xmax=lane_idx + 0.4, 
                              colors='#FFA54F', linewidth=10, alpha=0.2)
                    # 밴드 위에 수치 텍스트 표시
                    ax.text(lane_idx, bp * 1.05, f"{bp}bp", color='#FF7F24', 
                            fontsize=9, ha='center', fontweight='bold', va='bottom')

            # --- 결과 출력 ---
            st.pyplot(fig)

            # 다운로드 버튼
            buf = io.BytesIO()
            fig.savefig(buf, format="png", facecolor='#361F00')
            st.download_button(
                label="📥 Download Realistic Gel Image", 
                data=buf.getvalue(), 
                file_name="realistic_gel_output.png", 
                mime="image/png"
            )
    else:
        st.warning("위 입력창에 bp 값을 입력해 주세요 (예: 500, 1200).")
with tab4:
    st.header("Concentration Calculator")
    st.subheader("cDNA Synthesis Calculation")
    st.info("RNA 2,000ng(2μg)을 기준으로 cDNA 합성에 필요한 RNA와 Water의 양을 계산합니다.")

    # 1. 입력 섹션 (들여쓰기 없이 깔끔하게 배치)
    rna_conc = st.number_input(
        "현재 측정된 RNA 농도를 입력하세요 (ng/μL)", 
        min_value=0.0, 
        value=481.095, 
        format="%.3f",
        help="NanoDrop 등으로 측정된 값을 입력하세요."
    )

    # 2. 계산 로직
    if rna_conc > 0:
        # RNA 볼륨 = 2000 / RNA 농도
        rna_volume = 2000 / rna_conc
        
        # Water 볼륨 = 11 - RNA 볼륨
        water_volume = 11 - rna_volume

        # 3. 결과 출력
        st.markdown("---")
        res_col1, res_col2 = st.columns(2)
        
        with res_col1:
            st.metric(label="Required RNA Volume", value=f"{rna_volume:.2f} μL")
            
        with res_col2:
            if water_volume >= 0:
                st.metric(label="Nuclease-Free Water", value=f"{water_volume:.2f} μL")
            else:
                st.error("RNA 농도가 너무 낮아 11μL를 초과합니다!")
                st.warning("RNA 양을 줄이거나 샘플을 농축해야 합니다.")

        # 4. 프로토콜 가이드 (도움말)
        if water_volume >= 0:
            st.success(f"**Total Mixture (11μL) = RNA {rna_volume:.2f}μL + Water {water_volume:.2f}μL**")
            with st.expander("상세 믹스 구성 보기"):
                st.write(f"""
                1. RNA Sample: {rna_volume:.2f} μL (2,000ng 기준)
                2. N.F. Water: {water_volume:.2f} μL
                3. **Total (Template): 11.00 μL**
                ---
                *이후 RT Master Mix 등을 추가하여 최종 볼륨을 맞추세요.*
                """)
    else:
        st.warning("농도를 입력해 주세요.")


with tab5:
    st.header("Gene Visualization & 3D Structure")
    st.info("PDB 구조를 조회하거나, 아미노산 서열로부터 3D 구조를 예측합니다.")

    # 1. 시각화 함수 정의
    def render_3d_structure(pdb_data, format='pdb'):
        view = py3Dmol.view(width=800, height=500)
        if format == 'pdb':
            view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.zoomTo()
        showmol(view, height=500, width=800)

    # 2. 입력 방식 선택 (라디오 버튼)
    input_mode = st.radio(
        "분석 방식 선택",
        ["PDB ID로 조회 (이미 알려진 구조)", "단백질 서열로 예측 (새로운 구조)"],
        horizontal=True
    )

    st.markdown("---")

    # --- 모드 1: PDB ID로 조회 ---
    if input_mode == "PDB ID로 조회 (이미 알려진 구조)":
        col1, col2 = st.columns([1, 2])
        with col1:
            target_pdb = st.text_input("PDB ID 입력", value="4W5N", help="예: 4W5N(Argonaute), 1TUP(p53)").upper()
            search_btn = st.button("구조 불러오기", use_container_width=True)
        
        if search_btn or target_pdb:
            try:
                url = f"https://files.rcsb.org/view/{target_pdb}.pdb"
                response = requests.get(url)
                if response.status_code == 200:
                    st.success(f"PDB 데이터 로드 완료: {target_pdb}")
                    render_3d_structure(response.text)
                else:
                    st.error("해당 PDB ID를 찾을 수 없습니다.")
            except Exception as e:
                st.error(f"오류 발생: {e}")

    # --- 모드 2: 서열로 예측 (ESMFold API) ---
    else:
        st.subheader("AI Protein Structure Prediction")
        target_seq = st.text_area(
            "아미노산 서열(Amino Acid) 입력", 
            placeholder="예: MAG...", 
            height=150,
            help="Meta AI의 ESMFold를 사용하여 3D 구조를 예측합니다."
        )
        
        predict_btn = st.button("3D 구조 예측 시작 (AI Folding)", use_container_width=True)

        if predict_btn:
            if not target_seq or len(target_seq) < 5:
                st.warning("분석할 아미노산 서열을 입력해 주세요.")
            else:
                with st.spinner("AI가 단백질 구조를 접고 있습니다(Folding)... 약 10~30초 소요됩니다."):
                    try:
                        # Meta ESMFold API 사용
                        esmfold_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
                        response = requests.post(esmfold_url, data=target_seq.strip())
                        
                        if response.status_code == 200:
                            pdb_predicted = response.text
                            st.success("예측 완료!")
                            render_3d_structure(pdb_predicted)
                            
                            # 예측된 파일 다운로드 제공
                            st.download_button(
                                label="예측된 PDB 파일 다운로드",
                                data=pdb_predicted,
                                file_name="predicted_pwn_protein.pdb",
                                mime="text/plain"
                            )
                        else:
                            st.error("예측 실패. 서열이 너무 길거나 서버 응답에 문제가 있습니다.")
                    except Exception as e:
                        st.error(f"예측 중 오류 발생: {e}")


with tab6:
    st.header("Gene 조정")
    st.write("전사 및 번역을 실시")
    st.markdown("---") 

    # 서열 입력창 (들여쓰기 없이 배치)
    input_seq = st.text_area("분석할 DNA 서열을 입력하세요 (ATGC...)", height=200, placeholder="여기에 서열을 붙여넣으세요.", key="tab6_input")

    if input_seq:
        # 비정상 문자 제거 및 대문자 변환
        clean_seq = "".join(input_seq.split()).upper()
        my_seq = Seq(clean_seq)

        st.subheader("변환 결과")
        
        # 3단 컬럼 배치 (전사, 번역, 역상보)
        res_c1, res_c2, res_c3 = st.columns(3)
        
        with res_c1:
            st.info("**Transcription (DNA → RNA)**")
            rna_seq = my_seq.transcribe()
            st.code(str(rna_seq), language="text")
            st.download_button("RNA 다운로드", str(rna_seq), file_name="transcription.txt")

        with res_c2:
            st.success("**Translation (DNA → Protein)**")
            try:
                # 중간 정지 코돈(*)이 나올 수 있으므로 처리
                prot_seq = my_seq.translate()
                st.code(str(prot_seq), language="text")
                st.download_button("단백질 다운로드", str(prot_seq), file_name="translation.txt")
            except Exception as e:
                st.error("번역 실패: 서열 길이를 확인하세요 (3의 배수).")

        with res_c3:
            st.warning("🔄 **Reverse Complement**")
            rev_seq = my_seq.reverse_complement()
            st.code(str(rev_seq), language="text")
            st.download_button("역상보 서열 다운로드", str(rev_seq), file_name="rev_complement.txt")

     





