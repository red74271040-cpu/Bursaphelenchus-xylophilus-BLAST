#모듈 불러오기
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

      
         
with tab1:
    st.header("B. xylophilus 유전자 정보 검색")
    st.info("구축된 로컬 CDS 데이터베이스를 활용하여 재선충 유전자를 검색합니다.")

    # 좌우 레이아웃 분할 (입력창과 설명창)
    c1_in, c1_gui = st.columns([3, 2])
    
    with c1_in:
        st.subheader("Sequence Search")
        # 검색할 서열 입력 (프라이머 또는 유전자 조각)
        query_seq = st.text_area("dsRNA Primer Sequence input", height=150, 
                                 placeholder="분석할 서열(ATGC...)을 입력하십시오.", 
                                 key="pwn_local_search")
        
        if st.button("RUN LOCAL BLAST", use_container_width=True):
            if not query_seq or len(query_seq) < 15:
                st.warning("15bp 이상의 서열을 입력해 주십시오.")
            else:
                base_path = os.getcwd() 
                temp_query = os.path.join(base_path, "temp_query.fa")
                with open(temp_query, "w") as f:
                    f.write(f">Query\n{query_seq}")
                db_folder = os.path.join(base_path, "pwn_db")
                db_path = os.path.join(db_folder, "pwn_db")
                result_csv = os.path.join(base_path, "blast_result.csv")
        
        
                
                with st.spinner("로컬 데이터베이스 검색 중..."):
                    try:
                        # 3. BLAST 실행
                        # 만약 환경 변수 등록이 안 되어 있다면 blast_exe 경로를 직접 넣습니다.
                        cmd = [
                            "blastn", # 또는 blast_exe (경로를 직접 넣으려면)
                            "-query", temp_query,
                            "-db", db_path,
                            "-out", result_csv,
                            "-outfmt", "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                            "-task", "blastn-short"
                        ]
                        
                        subprocess.run(cmd, check=True, shell=True)
                        
                        # 4. 파일 크기 체크 (os.path.getsize로 수정)
                        if os.path.exists(result_csv) and os.path.getsize(result_csv) > 0:
                            df = pd.read_csv(result_csv, names=[
                                "Query", "Locus ID", "Identity(%)", "Length", 
                                "Mismatch", "Gaps", "Q_Start", "Q_End", 
                                "S_Start", "S_End", "E-value", "BitScore"
                            ])
                            st.success(f" 검색 완료: {len(df)}개의 결과를 찾았습니다.")
                            st.dataframe(df, use_container_width=True)
                        else:
                            st.error(" 매칭되는 결과가 없습니다. 서열을 확인하세요.")
                            
                    except Exception as e:
                        st.error(f"실행 중 오류 발생: {e}")

    with c1_gui:
        st.subheader("도우미")
        st.write("이 시스템은 **C:\\Users\\SAMSUNG\\blast_work** 폴더의 로컬 인덱스를 사용합니다.")
        with st.expander("결과 항목 설명"):
            st.markdown("""
            - **Locus ID**: 재선충 유전자 식별 번호 (예: BXYJ_LOCUS...)
            - **Identity**: 서열 일치율 (%)
            - **E-value**: 0에 가까울수록 통계적으로 유의미한 결과
            - **S_Start/End**: 재선충 유전자 내 매칭 위치
            """)
    # 기존 도우미 내용 하단에 추가
    st.markdown("---")
    st.subheader("🌐 NCBI Protein Quick Search")
    st.write("검색된 **Locus ID**를 입력하면 NCBI Protein DB로 바로 연결됩니다.")

    # 사용자로부터 ID 입력받기
    search_id = st.text_input("Enter Locus ID (e.g., BAE48369.1)", placeholder="Locus ID 입력...")

    if search_id:
        # NCBI Protein 검색 URL 생성
        ncbi_url = f"https://www.ncbi.nlm.nih.gov/protein/{search_id}"
        
        # 버튼 클릭 시 해당 페이지 열기 (Markdown 활용)
        st.markdown(f"""
            <a href="{ncbi_url}" target="_blank">
                <button style="
                    width: 100%;
                    background-color: #1a2a6c;
                    color: white;
                    padding: 10px;
                    border: none;
                    border-radius: 5px;
                    cursor: pointer;
                    font-weight: bold;
                ">NCBI에서 {search_id} 상세 정보 확인하기 ↗️</button>
            </a>
            """, unsafe_allow_html=True)
        
        st.caption("※ 새 창에서 NCBI Protein 데이터베이스가 열립니다.")

with tab2:
    st.header("si-Fi RNAi 분석 엔진")
    st.info("모드 선택에 따라 siRNA 효율성 측정 또는 타 생물군 오프타겟 위험도를 분석합니다.")

    # 1. 모드 및 설정 선택
    col_mode, col_set = st.columns([1, 1])
    with col_mode:
        mode_selection = st.radio("분석 모드 선택", ["RNAi Design Mode", "Off-target Prediction Mode"])
        mode = 0 if mode_selection == "RNAi Design Mode" else 1
    
    with col_set:
        st.write("⚙️ **Pipeline Parameters**")
        si_size = st.number_input("siRNA Size", value=21)
        mismatch_limit = st.slider("Allowed Mismatches", 0, 3, 0)

    # 2. 데이터베이스 및 서열 입력
    st.markdown("---")
    c1, c2 = st.columns(2)
    with c1:
        # 디자인 모드면 표적 유전자, 오프타겟 모드면 대조군 DB 선택 용도
        db_name = st.text_input("Database Name (in blast_work folder)", value="pwn_db")
        target_fasta = st.text_area("Target mRNA Sequence (FASTA)", height=200, placeholder=">Gene_ID\nATGC...")
    
    with c2:
        st.write("🔬 **Efficiency Filters**")
        use_strand = st.checkbox("Strand Selection (Energy)", value=True)
        use_end = st.checkbox("End Stability", value=True)
        use_access = st.checkbox("Target Site Accessibility", value=True)
        acc_window = st.number_input("Accessibility Window", value=8)

    # 3. 파이프라인 실행 버튼
    if st.button("RUN si-Fi PIPELINE", use_container_width=True):
        if not target_fasta.strip():
            st.error("분석할 서열을 입력해 주세요.")
        else:
            # 워킹 디렉토리 설정 (사용자 경로)
            work_dir = r"C:\Users\SAMSUNG\blast_work"
            if not os.path.exists(work_dir):
                os.makedirs(work_dir)

            with st.spinner("파이프라인 가동 중... (Bowtie & RNAplfold 실행)"):
                try:
                    # (1) 입력 서열 임시 저장
                    tmp_query = os.path.join(work_dir, "query_tmp.fa")
                    with open(tmp_query, "w") as f:
                        f.write(target_fasta)

                    # (2) SifiPipeline 클래스 호출
                    # sifi_pipeline.py에 정의된 클래스를 가져옵니다.
                    from sifi_pipeline import SifiPipeline
                    
                    pipeline = SifiPipeline(
                        bowtie_db=db_name,
                        db_location=work_dir + "\\",
                        query_sequences=tmp_query,
                        sirna_size=si_size,
                        mismatches=mismatch_limit,
                        accessibility_check=use_access,
                        accessibility_window=acc_window,
                        rnaplfold_location=work_dir,
                        bowtie_location=work_dir,
                        mode=mode, # 0: Design, 1: Off-target
                        strand_check=use_strand,
                        end_check=use_end,
                        end_stability_treshold=1.0,
                        target_site_accessibility_treshold=0.1,
                        temp_location=work_dir,
                        terminal_check=True,
                        no_efficience=False
                    )

                    # (3) 결과 도출
                    results = pipeline.run_pipeline
                    
                    if results:
                        img, json_path, table_data, main_dict, off_dict, err = results
                        
                        if err:
                            st.error(f"분석 중 오류: {err}")
                        else:
                            st.success("✅ 분석 완료!")
                            
                            # 디자인 모드: 효율성 중심 결과 출력
                            if mode == 0:
                                st.subheader("🎯 RNAi Design Results")
                                # general_helpers를 이용해 요약 테이블 생성
                                summary = general_helpers.get_table_data(json_path)
                                df_sum = pd.DataFrame(summary, columns=["Target", "Total hits", "Efficient hits"])
                                st.dataframe(df_sum, use_container_width=True)
                            
                            # 오프타겟 모드: 위험군 중심 결과 출력
                            else:
                                st.subheader("⚠️ Off-target Risk Analysis")
                                summary = general_helpers.get_table_data(json_path)
                                df_off = pd.DataFrame(summary, columns=["Organism/Gene", "Match Count", "High-Risk Hits"])
                                st.dataframe(df_off, use_container_width=True)

                            # 시각화 이미지 (show_plot.py 결과)
                            if img and os.path.exists(img):
                                st.image(img, caption="si-Fi Result Visualization")

                            # 상세 데이터 다운로드 (JSON 기반)
                            raw_data = general_helpers.prepare_json_data(json_path)
                            st.download_button("📥 상세 데이터(JSON) 저장", json.dumps(raw_data), "result.json")

                except Exception as e:
                    st.error(f"시스템 에러: {e}")
                    st.info("참고: bowtie.exe와 RNAplfold.exe가 C:\\Users\\SAMSUNG\\blast_work 폴더에 있어야 합니다.")
 
    




import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter, NullFormatter

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

     





